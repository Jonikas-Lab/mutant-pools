#!/usr/bin/env python
""" Join the by-mutant counts from multiple files generated by mutant_count_alignments.py into a single tab-separated file, along with position and gene information for each mutant.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: mutant_join_datasets.py [options] infile1 [infile2 infile3 ...] outfile """

# basic libraries
import sys, os, time
# other libraries
from collections import defaultdict
import unittest
# my modules
from general_utilities import write_header_data
import mutant_analysis_classes
from testing_utilities import run_functional_tests


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    test_folder = "test_data"
    dataset1 = "test_data/INPUT_mutants1_no-genes.txt"
    dataset2 = "test_data/INPUT_mutants2_with-genes.txt"

    tests = [("join-datasets__basic", "-o position %s %s -q"%(dataset1, dataset2)), 
             ("join-datasets__with-names", "-D dataset1,dataset2 -o position %s %s -q"%(dataset1, dataset2)),
             ("join-datasets__other-order", "-D dataset2,dataset1 -o position %s %s -q"%(dataset2, dataset1)),
            ]
    # MAYBE-TODO add run test for -A/-a option? 
    # MAYBE-TODO add run-tests for -X, -z, -Z?

    parser = define_option_parser()
    argument_converter = lambda parser,options,args: (args[:-1], args[-1], options)
    return run_functional_tests(tests, parser, main, test_folder, 
                                argument_converter=argument_converter, append_to_outfilenames='.txt') 


class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__(self):
        print "NO UNIT-TESTS FOR THIS MODULE"


######### Main function code #########

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    ### functionality options
    parser.add_option('-D', '--dataset_names', default=None, metavar='A,B,C,...', 
                      help="Comma-separated list of short dataset names (no spaces) to use in header "
                          +"(if none, dataset names will be derived from filenames) (default %default)")

    parser.add_option('-o', '--sort_data_key', choices=['position','read_count','none'], default='position', 
                      metavar='position|read_count|none', help="Sort the output data: by alignment position, read count, "
                         +"or don't sort at all (default %default) - sorting may be slow for large datasets!")

    parser.add_option('-A', '--gene_annotation_file', default=None, metavar='FILE', 
                      help="Tab-separated file to use to look up gene names/descriptions from IDs (default %default)")
    parser.add_option('-a', '--annotation_file_is_standard', action='store_true', default=False,
                      help="Use if file provided in -A is the standard Cre type (and missing a header) (default %default)")

    parser.add_option('-X', '--remove_mutants_from_file', metavar='FILE',
                      help='Remove all mutants present in FILE from the datasets (see -Y for read count cutoff).')
    parser.add_option('-z', '--remove_mutants_readcount_min', type='int', default=1, metavar='M',
                      help='When applying -X, only remove mutants with at least N reads in FILE (default %default).')
    parser.add_option('-Z', '--remove_mutants_min_is_perfect', action='store_true', default=False,
                      help='When applying -X with -z M, compare M to perfect readcount, not total. (default %default).')

    # TODO add --dont_count_cassette and --dont_count_other options like in mutant_count_alignments.py? 

    ### MAYBE-TODO add options to specify which MUTANTS to include, based on:
    #   - whether they show up in any/all datasets etc
    #   - perfect/imperfect read %
    #parser.add_option('-m', '--output_only_shared_mutants', action='store_true', default=False,
    #                  help="Only output the mutants that have non-zero counts in ALL input files (default %default)")
    #parser.add_option('-M', '--output_all_mutants', action='store_false', dest='output_only_shared_mutants',
    #                  help="Output all mutants, including ones that only appear in one input file (turns -o off).")
    #parser.add_option('-P', '--mutants_perfect_percent', type='int', default=0, metavar='K',
    #                  help="Include only mutants with K% or more perfectly aligned reads in all files (default %default)")
    #parser.add_option('-I', '--mutants_imperfect_percent', type='int', default=0, metavar='L',
    #                  help="Include only mutants with L% or more imperfectly aligned reads"
    #                      +"(calculated out of readcount sum over all infiles) (default %default)")
    # TODO move the -P/-I functionality to mutant_make_plots.py!
    #  Note: this -P/-I thing isn't as simple as it looks! 
    #   1) I'm not sure how to implement it - check the K% or L% for all files, any file, sum over all files?  Make sure it's consistent somehow, and that one mutant will always be in exactly one set between -P and -I (so either both use the sum over all files, or one requires its condition for all files and one for any file)  
    # OR just do it per file here, allowing some mutants to be present in one file but 0 in another, and let them get plotted like that? MIGHT BE THE BEST WAY.  
    # ANYWAY, in the final version, I think the joint-mutant file should just contain both total and perfect reads (possibly with options to do otherwise), and the plotting program should get all these -a/-p/-i/-P/-I options.
    #   2) As this is set up right now, I can't do that anyway, because for each infile there's only one value saved, either total or perfect or imperfect readcount - no way of getting a ratio. I'd have to rewrite a fair amount to fix that, and if I'm rewriting already, I should just do a full rewrite to use the new multi-dataset Insertional_mutant option from mutant_analysis_classes.py!

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


# MAYBE-TODO why is this so slow on large datasets?  (is that still the case after the rewrite?)
def main(infiles, outfile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing. 
    Print final dataset to outfile (if given); return final multi-dataset object and the list of dataset names in order.
    The options argument should be generated by an optparse parser.
    """

    # parse all infiles, print summaries to stdout if requested
    all_datasets = {}

    if options.dataset_names:
        dataset_names = options.dataset_names.split(',')
        if not len(dataset_names)==len(infiles):
            raise ValueError("If dataset names are provided via -D option, you must provide the same number of names "
                             +"as the total number of infiles! We have %s names and %s infiles."%(len(dataset_names), 
                                                                                                  len(infiles)))
    else:
        dataset_names = [os.path.splitext(os.path.basename(infile))[0] for infile in infiles]

    for dataset_name,infile in zip(dataset_names,infiles):
        if options.verbosity_level>1:   print "parsing input file %s - time %s."%(infile, time.ctime())
        current_dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset(infile=infile)
        current_dataset.count_adjacent_mutants(OUTPUT=None)
        # TODO once read_data_from_file actually reads mutant-merging info, get the adjacent-max-distance from that and pass that to current_dataset.count_adjacent_mutants instead of using the default value
        all_datasets[dataset_name] = current_dataset
        if options.verbosity_level>0:   print "%s mutants in dataset from input file %s"%(len(current_dataset), infile)
        elif options.verbosity_level>1: current_dataset.print_summary()
    
    # merge datasets into one multi-dataset object
    if options.verbosity_level>1:   print "merging the mutant data into combined dataset - time %s."%(time.ctime())
    multi_dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset(multi_dataset=True)
    multi_dataset.populate_multi_dataset(all_datasets, overwrite=False, check_gene_data=True)
    # make sure the datasets are in the same order as they were given on the command-line
    #  (using all_datasets to initialize multi_dataset didn't give an order, since all_datasets is a dictionary)
    multi_dataset.dataset_order = dataset_names
    # print varying amounts of summary data to stdout
    if options.verbosity_level>0:   print "total %s mutants present in combined dataset"%(len(multi_dataset))
    elif options.verbosity_level>0: multi_dataset.print_summary()

    ### optionally remove mutants based on another dataset
    if options.remove_mutants_from_file:
        other_dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset(infile=options.remove_mutants_from_file)
        multi_dataset.remove_mutants_based_on_other_dataset(other_dataset, 
                 readcount_min=options.remove_mutants_readcount_min, perfect_reads=options.remove_mutants_min_is_perfect)

    # if requested, add gene annotation info from separate file
    if options.gene_annotation_file:
        if options.verbosity_level>1: 
            print "adding gene annotation from file %s - time %s."%(options.gene_annotation_file, time.ctime())
        multi_dataset.add_gene_annotation(options.gene_annotation_file, 
                                               if_standard_Cre_file=options.annotation_file_is_standard)

    # print full data to outfile, unless there is no outfile name given
    if outfile:
        if options.verbosity_level>1:   
            print "printing combined dataset output to file %s - time %s."%(outfile, time.ctime())
        with open(outfile,'w') as OUTFILE:
            write_header_data(OUTFILE,options)
            OUTFILE.write("### DATASET SUMMARIES:\n")
            multi_dataset.print_summary(OUTPUT=OUTFILE, line_prefix="#  ", header_prefix="## ")
            OUTFILE.write("### HEADER AND DATA:\n")
            multi_dataset.print_data(OUTPUT=OUTFILE, sort_data_by=options.sort_data_key, header_line=True)

    return multi_dataset, dataset_names


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # if ran with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        print("\n * unit-tests for the mutant_analysis_classes.py module")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(mutant_analysis_classes)
        unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        print("\n * unit-tests for this module (%s)"%sys.argv[0])
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise parse the arguments and run main function
    if len(args)<2:
        parser.print_help()
        sys.exit("\nError: at least one infile and exactly one outfile are required!")
    outfile = args[-1]
    infiles = args[:-1]

    main(infiles, outfile, options)

