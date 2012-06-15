#!/usr/bin/env python
""" Take a deepseq alignment file; group the reads into insertional mutants.
Output a line-per-mutant file containing position info, gene annotation data (optional), total/perfect read count, number of distinct sequences, and some of the sequences/counts (optional), and a summary of read/mutant counts etc. 
Also output a line-per-gene file containing gene ID/name/position/annotation and the number and read-counts of mutants that were inserted into that gene (optional, NOT IMPLEMENTED YET).
Output files are in simple tab-separated plaintext format, with all header/summary lines starting with #.

Grouping reads into mutants is currently done by alignment position (other options such as sequence clustering or grouping 1bp-adjacent positions together may be implemented later). 

The input file should be a SAM-format deepseq alignment file created by bowtie, novoalign, or other deepseq aligner programs (tested mainly on bowtie), with optionally a metadata file created by my deepseq_alignment_wrapper.py or deepseq_preprocessing_wrapper.py scripts.
The program assumes the SAM file contains only unique matches (i.e. each read was reported as aligning to at most one genomic location).
It also only deals with single-end alignments at the moment.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: mutant_count_alignments.py [options] infile1 [infile2 infile3 ...] outfile """

# basic libraries
import sys, os, time
import filecmp
import unittest
# other libraries
import HTSeq
# my modules
from general_utilities import write_header_data
from testing_utilities import run_functional_tests
import mutant_analysis_classes

def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    test_folder = "test_data"
    aln_infile1 = "test_data/INPUT_alignment1.sam"
    aln_infile2 = "test_data/INPUT_alignment2_for-genes.sam"
    aln_infile3 = "test_data/INPUT_alignment3_for-merging.sam"
    gff_genefile = "test_data/INPUT_gene-data.gff3"
    dataset_to_remove = "test_data/count-aln__merge-adjacent2-r3.txt"

    test_runs = [('cassette-end-5prime', "-e 5prime -r forward -U -n3", [aln_infile1]),
                 ('cassette-end-3prime', "-e 3prime -r forward -U -n3", [aln_infile1]),
                 ('read-direction-reverse', "-r reverse -e 5prime -U -n3", [aln_infile1]),
                 ('u_unknown-not-as-match', "-u -e 5prime -r forward -n3", [aln_infile1]),
                 ('dont-count-cassette', "-I -e 5prime -r forward -U -n3", [aln_infile1]),
                 ('ignore-cassette', "-i -e 5prime -r forward -U -n3", [aln_infile1]),
                 ('sorted-by-count', "-o read_count -e 5prime -r forward -U -n3", [aln_infile1]),
                 ('with-gene-info', "-e 5prime -r forward -U -g %s -d -n0"%gff_genefile, [aln_infile2]),
                 ('multiple-infiles', "-e 5prime -r forward -U -n0", [aln_infile1,aln_infile2]),
                 ('merge-adjacent-none', "-n0", [aln_infile3]),
                 ('merge-adjacent1-r3', "-M --merge_max_distance=1 --merge_count_ratio=3 -n0", [aln_infile3]),
                 ('merge-adjacent1-r1', "-M --merge_max_distance=1 --merge_count_ratio=1 -n0", [aln_infile3]),
                 ('merge-adjacent2-r3', "-M --merge_max_distance=2 --merge_count_ratio=3 -n0", [aln_infile3]), 
                 ('remove-from-other-all', "-X %s -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-from-other-min4', "-X %s -z4 -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-from-other-perfect', "-X %s -Z -z4 -n0"%dataset_to_remove, [aln_infile2]),
                ]
    # MAYBE-TODO add run-test for --gene_annotation_file?
    # MAYBE-TODO add run-test for --input_collapsed_to_unique?  Or is that unit-tested already?


    # convert tests into (testname, arg_and_infile_string) format, adding the options that are always used
    test_names_and_args = [('count-aln__'+testname, test_args+' -q '+' '.join(infiles)) 
                           for testname,test_args,infiles in test_runs]

    parser = define_option_parser()
    argument_converter = lambda parser,options,args: (args[:-1], args[-1], options)
    return run_functional_tests(test_names_and_args, parser, main, test_folder, 
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

    # taken:     aAbBcCdDe---g-h-iI------mMn-o---q-r---tTuUvVwWxX-YzZ  
    # free:      ---------EfF-G-H--jJkKlL---N-OpP-Q-RsS----------y---  

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    ### functionality options
    parser.add_option('-e', '--read_cassette_end', choices=mutant_analysis_classes.SEQ_ENDS, default='5prime', 
                      metavar='|'.join(mutant_analysis_classes.SEQ_ENDS), 
                      help="Which end of the cassette are the sequenced reads from? (default %default).")
    parser.add_option('-r','--read_direction', choices=mutant_analysis_classes.SEQ_DIRECTIONS, default='forward',
                      metavar='|'.join(mutant_analysis_classes.SEQ_DIRECTIONS), 
                      help="Is the read in the forward or reverse direction compared to the cassette? (default %default).")

    parser.add_option('-M', '--merge_adjacent_mutants', action="store_true", default=False, 
                      help="Merge adjacent mutants if they satisfy the -w/-W constraints (default %default)")
    parser.add_option('-w', '--merge_max_distance', type='int', default=1, metavar='N',
                      help="For -M: merge mutants only if they're at most N bases distant (default %default)")
    parser.add_option('-W', '--merge_count_ratio', type='int', default=1000, metavar='K',
                      help="For -M: merge mutants only if one has K x fewer reads than the other (default %default)")

    parser.add_option('-X', '--remove_mutants_from_file', metavar='FILE',
                      help='Remove all mutants present in FILE from the datasets (see -z/-Z for read count cutoff).')
    parser.add_option('-z', '--remove_mutants_readcount_min', type='int', default=1, metavar='M',
                      help='When applying -X, only remove mutants with at least N reads in FILE (default %default).')
    parser.add_option('-Z', '--remove_mutants_min_is_perfect', action='store_true', default=False,
                      help='When applying -X with -z M, compare M to perfect readcount, not total. (default %default).')

    # extremely minor functionality options, do we even care??
    parser.add_option('-u', '--treat_unknown_as_match', action="store_true", default=False, 
                      help="When counting perfect reads, treat undefined alignment regions as matches (default %default)")
    parser.add_option('-U', '--dont_treat_unknown_as_match', action="store_false", dest='treat_unknown_as_match',
                      help="Turn -u off.")
    # MAYBE-TODO add user-provided mutation-count cutoffs like in old deepseq_count_alignments.py, instead of just all reads and perfet reads?   Currently useless, since we're only allowing one mutation in bowtie.  parser.add_option('-m', '--mutation_cutoffs', default="1,3,10", metavar="<comma-separated-int-list>")

    ### input options
    parser.add_option('-c','--input_collapsed_to_unique', action='store_true', default=False, 
                      help="Use to get correct original total read counts if the data was collapsed to unique sequences "
                          +"using fastx_collapser before alignment (default %default).")
    parser.add_option('-C','--input_not_collapsed_to_unique', action='store_false', dest="input_collapsed_to_unique", 
                      help="Turn -c off.")
    parser.add_option('-m', '--input_metadata_file', default='AUTO', metavar='FILE', 
                      help="File containing preprocessing and alignment metadata (scripts/options used etc). "
                          +"Can be a filename, AUTO for <infile_basename>_info.txt (warning will be raised if not found), "
                          +"or NONE to not look for a metadata file at all. Default %default.")

    ### gene-finding and gene-annotation options 
    parser.add_option('-g', '--gene_position_reference_file', default=None, metavar='FILE', 
                      help="File to use to look up gene IDs based on chromosomal location (default %default)")
    parser.add_option('-d', '--detailed_gene_features', action="store_true", default=True,
                      help="Find out what part of the gene (UTR,intron,exon) a mutant hit, based on the -g file "
                          +"(default %default). May take a lot of memory - increase -Y option value to fix that.")
    parser.add_option('-D', '--no_detailed_gene_features', action="store_false", dest='detailed_gene_features',
                      help="Turns -d off.")
    parser.add_option('-Y', '--N_detail_run_groups', type="int", default=5, metavar='N', 
                      help="How many passes to split reading the detailed_gene_features into (default %default) "
                          +"- may take a lot of memory (and CPU) if read in a single pass; too many passes waste CPU.")
    # MAYBE-TODO add a "flank" option (with variable size), to catch mutants that are in the flanks of genes? Do we care?
    # MAYBE-TODO add a "negative flank" option (with variable size), to ignore mutants that are in the start/end of genes?

    parser.add_option('-A', '--gene_annotation_file', default=None, metavar='FILE', 
                      help="Tab-separated file to use to look up gene names/descriptions from IDs (default %default)")
    parser.add_option('-a', '--annotation_file_is_standard', action='store_true', default=False,
                      help="Use if file provided in -A is the standard Cre type and (missing a header) (default %default)")

    ### output format options
    parser.add_option('-n', '--N_sequences_per_group', type='int', default=2, metavar='N', 
                      help="How many most common sequences should be shown per group? (default %default)")
    parser.add_option('-o', '--sort_data_key', choices=['position','read_count','none'], default='position', 
                      metavar='position|read_count|none', help="Sort the output data: by alignment position, read count, "
                         +"or don't sort at all (default %default) - sorting may be slow for large datasets!")
    parser.add_option('-I', '--dont_count_cassette', action='store_true', default=False, 
                      help="Count cassette reads and print the count in the header; "
                          +"otherwise treat them normally. (default %default) (also see -i)")
    parser.add_option('-i', '--ignore_cassette', action='store_true', default=False,
                      help="Ignore reads aligning to cassette (just print total count in the header as removed) "
                          +"(default %default) (also see -E)")

    parser.add_option('-V', '--verbosity_level', action="store_true", default=1, 
                      help="How much information to print to STDOUT: 0 - nothing, 1 - summary only, "
                          +"2 - summary and progress reports. (Default %default).")
    parser.add_option('-q', '--quiet', action="store_const", const=0, dest='verbosity_level', help="Equivalent to -V 0.")
    parser.add_option('-v', '--verbose', action="store_const", const=2, dest='verbosity_level', help="Equivalent to -V 2.")

    return parser



def add_discarded_reads_from_metadata_file(infiles, input_metadata_file, verbosity_level):
    """ Parse metadata files to get total discarded read count; return None if it cannot be determined. 
    """
    # if the option specified no metadata files, total discarded readcount cannot be determined
    if input_metadata_file == 'NONE':
        return None
    # make sure the -m option has a value that will work with the number of infiles
    if len(infiles)>1 and not input_metadata_file=='AUTO':
        print "Warning: when multiple input files are given, the -m option must be NONE or AUTO - ignoring other value."
        return None
    # get the discarded read count for each infile; only return the total at the end, if all values are found
    discarded_counts = []
    for infile in infiles:
        # take the actual input_metadata_file value as the file name, or infer it from the infile name if 'AUTO'
        if input_metadata_file == 'AUTO':
            curr_input_metadata_file = os.path.splitext(infile)[0] + '_info.txt'
            if verbosity_level>1:  
                print 'Automatically determining metadata input file name: %s'%curr_input_metadata_file
        else:
            curr_input_metadata_file = input_metadata_file
            if verbosity_level>1:  
                print 'Metadata input file name provided in options: %s'%curr_input_metadata_file
        # if file is missing, total discarded readcount cannot be determined
        if not os.path.exists(curr_input_metadata_file):
            if verbosity_level>0:
                print 'Warning: metadata input file %s not found! Proceeding without it.'%curr_input_metadata_file
            return None
        # go through the metadata file to find the line with the discarded read count; 
        #  if the line isn't found, total discarded readcount cannot be determined
        line_found = False
        for line in open(curr_input_metadata_file):
            if line.startswith('## reads removed: '):
                discarded_counts.append(int(line.split()[3]))
                line_found = True
                break
        if not line_found:
            if verbosity_level>0:
                print("Warning: metadata input file %s didn't contain discarded read count line! "%curr_input_metadata_file
                      +"Proceeding without it.")
            return None
    # if some metadata files were missing the discarded counts line, total discarded readcount cannot be determined
    if len(discarded_counts)<len(infiles):
        if verbosity_level>0:
            print "Warning: discarded read count not found for some files! Ignoring all values, keeping 'unknown'."
        return None
    # if the discarded read counts for all infiles were found, return the sum as the total discarded read count
    return sum(discarded_counts)
    
    

def main(infiles, outfile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """
    ### parse/process/reformat some options
    options.count_cassette = not options.dont_count_cassette

    ### generate empty alignment set object with basic read position/orientation properties defined by options
    all_alignment_data = mutant_analysis_classes.Insertional_mutant_pool_dataset(options.read_cassette_end, 
                                                                                 options.read_direction=='reverse')

    ### parse preprocessing/alignment metadata file to get discarded read count, pass it to all_alignment_data
    #   (all_alignment_data initializes it to 'unkown', so if file is not given or can't be found, no need to do anything)
    N_discarded = add_discarded_reads_from_metadata_file(infiles, options.input_metadata_file, options.verbosity_level)
    if N_discarded is not None:     all_alignment_data.add_discarded_reads(N_discarded)
    # MAYBE-TODO also get the final total number of reads from the metadata infile and make sure it's the same 
    #   as the number of processed reads I get from all_alignment_data.print_summary()?
    
    ### parse input file and store data - the add_alignment_reader_to_data function here does pretty much all the work!
    for infile in infiles:
        # initialize a parser for the SAM infile
        if options.verbosity_level>1: print "parsing input file %s - time %s."%(infile, time.ctime())
        infile_reader = HTSeq.SAM_Reader(infile)
        # fill the new alignment set object with data from the infile parser
        all_alignment_data.add_alignment_reader_to_data(infile_reader, 
                                                        uncollapse_read_counts = options.input_collapsed_to_unique, 
                                                        treat_unknown_as_match = options.treat_unknown_as_match, 
                                                        count_cassette = options.count_cassette, 
                                                        ignore_cassette = options.ignore_cassette)

    ### optionally merge adjacent mutants (since they're probably just artifacts of indels during deepseq)
    if options.merge_adjacent_mutants: 
        all_alignment_data.merge_adjacent_mutants(merge_max_distance = options.merge_max_distance, 
                                                  merge_count_ratio = options.merge_count_ratio, 
                                                  dont_change_positions = True, 
                                                  verbosity_level=options.verbosity_level)

    ### optionally remove mutants based on another dataset
    if options.remove_mutants_from_file:
        other_dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset(infile=options.remove_mutants_from_file)
        all_alignment_data.remove_mutants_based_on_other_dataset(other_dataset, 
                 readcount_min=options.remove_mutants_readcount_min, perfect_reads=options.remove_mutants_min_is_perfect)

    ### optionally parse gene position/info files and look up the genes for each mutant in the data
    if options.gene_position_reference_file is not None:
        genefile = options.gene_position_reference_file
        if options.verbosity_level>1: print "adding genes from file %s to mutant data - time %s."%(genefile, time.ctime())
        all_alignment_data.find_genes_for_mutants(genefile, detailed_features=options.detailed_gene_features, 
                                                  ignore_cassette=options.ignore_cassette, 
                                                  N_run_groups=options.N_detail_run_groups, 
                                                  verbosity_level=options.verbosity_level)

        # if we have gene info, optionally also add annotation
        if options.gene_annotation_file:
            if options.verbosity_level>1: 
                print "adding gene annotation from file %s - time %s."%(options.gene_annotation_file, time.ctime())
            all_alignment_data.add_gene_annotation(options.gene_annotation_file, 
                                                   if_standard_Cre_file=options.annotation_file_is_standard)

    ### output
    # print summary info to stdout
    if options.verbosity_level>1:   print "\nDATA SUMMARY:"
    if options.verbosity_level>0:   all_alignment_data.print_summary()
    # print full data to outfile
    if options.verbosity_level>1:   print "printing output - time %s."%time.ctime()
    with open(outfile,'w') as OUTFILE:
        write_header_data(OUTFILE,options)
        OUTFILE.write("### SUMMARY:\n")
        all_alignment_data.print_summary(OUTFILE, line_prefix="#  ", header_prefix="## ")
        OUTFILE.write("### HEADER AND DATA:\n")
        all_alignment_data.print_data(OUTPUT=OUTFILE, sort_data_by=options.sort_data_key, 
                                      N_sequences=options.N_sequences_per_group, 
                                      header_line=True, header_prefix='# ')


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

