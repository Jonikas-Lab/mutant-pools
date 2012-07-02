#! /usr/bin/env python

"""
______
 -- Weronika Patena, 2012
USAGE: mutant_growth_rates.py [options] T0_infile T1_infile outfile
"""
# TODO write proper docstring!

# standard library
from __future__ import division
import sys
import unittest
from collections import defaultdict
# other packages
from numpy import log2, mean, median, std, isnan, isinf
# my modules
from parse_annotation_file import parse_gene_annotation_file
import mutant_join_datasets


def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    ### functionality options
    parser.add_option('-G', '--generations', type='int', default=7, metavar='N', 
                      help="The number of generations between the start and end timepoints (default %default).")
    parser.add_option('-R', '--pool_growth_rate', type='float', metavar='X', 
                      help="Overall growth rate of the culture (REQUIRED).")
    parser.add_option('-Z', '--minimum_growth_rate_zero', action='store_true', default=False,
                      help="Any growth rates below 0 will be set to 0 (default %default).")
    # TODO why is -Z and -m messed up weirdly when I run mutant_growth_rates.py -h on the command-line???

    parser.add_option('-m', '--readcount_min_T0', type='int', default=100, metavar='N', 
                      help="Minimum T0 readcount to calculate growthrate (default %default).")
    parser.add_option('-M', '--readcount_min_T1', type='int', default=0, metavar='N', 
                      help="Minimum T1 readcount to calculate growthrate (default %default).")

    parser.add_option('-A', '--gene_annotation_file', default=None, metavar='FILE', 
                      help="Tab-separated file to use to look up gene names/descriptions from IDs (default %default)")
    parser.add_option('-a', '--annotation_file_is_standard', action='store_true', default=False,
                      help="Use if file provided in -A is the standard Cre type (and missing a header) (default %default)")

    parser.add_option('-r', '--replace_missing_growthrates', type='float', default=float('nan'), metavar='X', 
                      help="What value to print when growthrate can't be calculated (default: %default).")
    parser.add_option('-D', '--dataset_names', default=None, metavar='A,B', 
                      help="Comma-separated list of short dataset names (no spaces) to use in header "
                          +"(if None, dataset names will be derived from filenames) (default %default)")

    ### command-line verbosity options
    parser.add_option('-V', '--verbosity_level', action='store_true', default=1, 
                      help="How much information to print to STDOUT: 0 - nothing, 1 - summary only, "
                          +"2 - summary and progress reports. (Default %default).")
    parser.add_option('-q', '--quiet', action='store_const', const=0, dest='verbosity_level', help="Equivalent to -V 0.")
    parser.add_option('-v', '--verbose', action='store_const', const=2, dest='verbosity_level', help="Equivalent to -V 2.")

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    return parser


def calculate_growth_rate(T0_readcount, T1_readcount, T0_total, T1_total, overall_growthrate, generations, 
                          T0_min, T1_min, minimum_value=None, NA_result=float('nan')):
    """ Calculate mutant growth rate: k = K * (1 + log2(x1/x0)/G); return NA_result if readcounts below T0_min and T1_min.
    x1 and x0 are the readcount of mutant x divided by the total sample readcount at timepoints T1 and T0,
    G is the number of doubling times the culture was grown for, and K is the overall growth rate of the culture.
    """
    if T0_readcount<T0_min or T1_readcount<T1_min:
        return NA_result
    T0_normalized = T0_readcount/T0_total
    T1_normalized = T1_readcount/T1_total
    growthrate = overall_growthrate * (1 + log2(T1_normalized/T0_normalized)/generations)
    if minimum_value is not None and growthrate<minimum_value:
        return minimum_value
    else:
        return growthrate


def main(args, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """

    try:
        infile_T0, infile_T1, outfile = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: exactly two infiles and exactly one outfile are required!")

    infiles = [infile_T0, infile_T1]

    # most of the work (making the multi-dataset with proper names, gene annotation, etc) is already implemented 
    #  in mutant_join_datasets.main - returns multi-dataset mutant_analysis_classes.Insertional_mutant_pool_dataset object.
    # use lower verbosity by one for the run, then reset it to old value.
    original_verbosity = options.verbosity_level
    options.verbosity_level = max(0, original_verbosity-1)
    options.remove_mutants_from_file = None
    multi_dataset, dataset_names = mutant_join_datasets.main(infiles, None, options)
    options.verbosity_level = original_verbosity
    T0_name, T1_name = dataset_names

    for mutant in multi_dataset:
        T0_readcount = mutant.by_dataset[T0_name].total_read_count
        T1_readcount = mutant.by_dataset[T1_name].total_read_count
        # TODO what should I be using as the total? processed_read_count? aligned_read_count? Something else?
        # TODO optionally use perfect instead of total read counts?
        T0_total = multi_dataset.summary[T0_name].processed_read_count
        T1_total = multi_dataset.summary[T1_name].processed_read_count
        mutant.growthrate = calculate_growth_rate(T0_readcount, T1_readcount, T0_total, T1_total, 
                                  options.pool_growth_rate, options.generations, 
                                  options.readcount_min_T0, options.readcount_min_T1, 
                                  0 if options.minimum_growth_rate_zero else None, options.replace_missing_growthrates)
        # TODO this is just the absolute growth rate - how about the relative one?  Add that too?
    # NOTE: If log2(x1/x2) is a negative value that's greater than G, we end up with a negative growth rate, because the equation doesn't take cell death into account.  Should we set 0 as the minimum?  Maybe... Worth looking at the raw data first, though.

    # TODO what's the median growth rate?  Some kind of histogram would be nice, really...
    # TODO I get -inf growthrate values - what to do with those?
    if options.verbosity_level>0:
        median_growthrate = median([m.growthrate for m in multi_dataset if not isnan(m.growthrate)])
        N_mutants = len([1 for m in multi_dataset if not isnan(m.growthrate)])
        print "median growthrate: %s (%s mutants)"%(median_growthrate, N_mutants)
        for cutoff in (0.2, 0.3, 0.5, 0.8):
            print "mutants with growthrate %s%% below and above median: %s, %s"%(100*cutoff, 
                                         len([m for m in multi_dataset if m.growthrate<=median_growthrate*(1-cutoff)]), 
                                         len([m for m in multi_dataset if m.growthrate>=median_growthrate*(1+cutoff)]))
            # TODO figure out if the number of mutants below and above the median is statistically significant?

    # print the data in some sensibly plottable way!
    with open(outfile, 'w') as OUTFILE:
        _print_growthrate_data_tmp(multi_dataset, T0_name, T1_name, OUTFILE)


def _print_growthrate_data_tmp(multi_dataset, T0_name, T1_name, OUTPUT=sys.stdout, header_prefix="# "):
    # modified copy of Insertional_mutant_pool_dataset.print_data - TODO replace this with something with less repetition!

    header = ['chromosome','strand','min_position','full_position', 'gene','orientation','feature','main_sequence']
    header += ['growth_rate', 'reads_T0','reads_T1']
    header += multi_dataset.gene_annotation_header
    OUTPUT.write(header_prefix + '\t'.join(header) + "\n")

    # create "empty" annotation line with the correct number of fields, for genes that weren't in annotation file
    if multi_dataset.gene_annotation_header:
        missing_gene_annotation_data = ['NO GENE DATA'] + ['' for x in multi_dataset.gene_annotation_header[:-1]]

    ### for each mutant, print the mutant data line (different for normal and multi-datasets)
    # need that if-else with an isnan treated as inf because nan's don't sort properly!
    for mutant in sorted(iter(multi_dataset), key = lambda m: float('inf') if isnan(m.growthrate) else m.growthrate):
        mutant_data = [mutant.position.chromosome, mutant.position.strand, mutant.position.min_position, 
                       mutant.position.full_position, mutant.gene, mutant.orientation, mutant.gene_feature] 
        mutant_data += [mutant.get_main_sequence(1)[0]]
        mutant_data += [mutant.growthrate, mutant.by_dataset[T0_name].total_read_count, 
                        mutant.by_dataset[T1_name].total_read_count]
        # add gene annotation, or a line with the right number of fields if gene annotation is missing
        if multi_dataset.gene_annotation_header:
            if mutant.gene_annotation:  mutant_data += mutant.gene_annotation
            else:                       mutant_data += missing_gene_annotation_data
        OUTPUT.write('\t'.join([str(x) for x in mutant_data]))
        OUTPUT.write('\n')


def _read_growthrate_data_from_file_tmp(self, infile, assume_new_sequences=False):
    # modified copy of Insertional_mutant_pool_dataset.read_data_from_file - TODO replace this with something with less repetition!
    if self.multi_dataset:  raise MutantError("read_data_from_file not implemented for multi-datasets!")
    for line in open(infile):
        # LATER-TODO get unaligned/discarded/etc read count from summary, so we can keep track of full counts!
        # ignore comment and header lines, parse other tab-separated lines into values
        if line.startswith('#'):                                        continue
        # this is needed only when reading test reference files
        if line.startswith('<REGEX>#') or line.startswith('<IGNORE>'):  continue       
        if line.startswith('chromosome\tstrand\tmin_position\t'):       continue
        fields = line.split('\t')
        chromosome = fields[0]
        strand = fields[1]
        min_pos = int(fields[2])
        full_pos = fields[3]
        gene, orientation, gene_feature = fields[4:7]
        total_reads,perfect_reads,sequence_variants = [int(x) for x in fields[7:10]]
        # generate new mutant if necessary; add counts and gene info to mutant (USE IMMUTABLE POSITIONS BY DEFAULT)
        position = Insertion_position(chromosome, strand, full_position=full_pos, immutable=True)
        curr_mutant = self.get_mutant(position)
        curr_mutant.add_counts(total_reads,perfect_reads,sequence_variants,assume_new_sequences)
        curr_mutant.update_gene_info(gene, orientation, gene_feature)
        # get however many specific sequences/counts are listed (this is variable)
        sequence_fields = fields[10::2]
        count_fields = fields[11::2]
        for seq, count in zip(sequence_fields, count_fields):
            if int(count)>0:
                assert seq!=''
                curr_mutant.sequences_and_counts[seq] += int(count)
        # add to dataset total read/mutant counts
        # MAYBE-TODO Might just want to write a function that recalculates all of the total counts below, to be ran at the end of read_data_from_file and add_alignment_reader_to_data and I guess find_genes_for_mutants, instead of doing it this way
        summ = self.summary
        summ.processed_read_count += total_reads
        summ.aligned_read_count += total_reads
        summ.perfect_read_count += perfect_reads
        summ.strand_read_counts[strand] += total_reads
        if gene==SPECIAL_GENE_CODES.not_found:        summ.mutants_not_in_genes += 1
        elif gene in SPECIAL_GENE_CODES.all_codes:    summ.mutants_undetermined += 1  # the two codes beside not_found
        else:                                         summ.mutants_in_genes += 1
        if orientation not in ['?','-']:              summ.mutant_counts_by_orientation[orientation] += 1
        if gene_feature not in ['?','-']:             summ.mutant_counts_by_feature[gene_feature] += 1
    # TODO update to deal with gene annotation fields?-
    


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    from testing_utilities import run_functional_tests
    test_folder = "test_data"
    sys.exit("NO TESTS DEFINED!")
    # tests in (testname, [test_description,] arg_and_infile_string) format
    test_runs = [ ]
    # argument_converter converts (parser,options,args) to the correct argument order for main
    argument_converter = lambda parser,options,args: (args, options)
    # use my custom function to run all the tests, auto-detect reference files, compare to output.
    return run_functional_tests(test_runs, define_option_parser(), main, test_folder, 
                                argument_converter=argument_converter, append_to_outfilenames='.txt') 
    # LATER-TODO add run-tests!


class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__(self):
        sys.exit("NO UNIT-TESTS FOR THIS MODULE")
    # LATER-TODO add unit-tests!


if __name__=='__main__':
    parser = define_option_parser()
    options,args = parser.parse_args()

    # if run with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        #print("\n * unit-tests for the ______ module")
        #test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(______
        #unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        print("\n * unit-tests for this module (%s)"%sys.argv[0])
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise pass the arguments to the main function
    main(args, options)
