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
import mutant_analysis_classes


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    test_infile1 = "test_data/test_input.sam"
    test_infile2 = "test_data/test_input2__for-genes.sam"
    test_constant_options = "-H 1 -s -n3 -o position -q"
    #  (using -s and -H 1 to get all relevant info but not have to deal with changing timestamps/etc)
    test_runs = [("-e 5prime -r forward -U", [test_infile1], "test_data/test_output__5prime.txt"),
                 ("-e 3prime -r forward -U", [test_infile1], "test_data/test_output__3prime.txt"),
                 ("-r reverse -e 5prime -U", [test_infile1], "test_data/test_output__reverse.txt"),
                 ("-u -e 5prime -r forward", [test_infile1], "test_data/test_output__u.txt"),
                 ("-b special_chromosome -e 5prime -r forward -U", [test_infile1], "test_data/test_output__b.txt"),
                 ("-B special_chromosome -e 5prime -r forward -U", [test_infile1], "test_data/test_output__B.txt"),
                 ("-o read_count -e 5prime -r forward -U", [test_infile1], "test_data/test_output__sorted-by-count.txt"),
                 ("-n 0 -e 5prime -r forward -U -g test_data/test_reference.gff3 -d", [test_infile2], 
                                                                      "test_data/test_output2__with-genes.txt"),
                 ("-n 0 -e 5prime -r forward -U", [test_infile1,test_infile2], "test_data/test_output12__multiple.txt")]
    #  Most common group: mutation_yes, position 204-?, + strand: 6 reads (20% of aligned)
    for variable_option_string, test_infiles, reference_file in test_runs:
        option_string = test_constant_options + ' ' + variable_option_string
        print(" * New test run, with options: %s (infiles %s; reference outfile %s)"%(option_string,
                                                                          ', '.join(test_infiles),reference_file))
        # regenerate options with test argument string
        parser = define_option_parser()
        (options, _) = parser.parse_args(option_string.split())
        outfile = "test_data/test_output.txt"
        run_main_function(test_infiles, outfile, options)
        # compare outfile to reference file: remove outfile and keep going if correct, otherwise exit with message.
        if filecmp.cmp(outfile, reference_file, shallow=False):
            os.remove(outfile)
        else:
            print("TEST FAILED!!  Reference file %s and output file %s differ - PLEASE COMPARE."%(reference_file,outfile))
            return 1
    print("*** Test runs finished - EVERYTHING IS FINE. ***")
    return 0


class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__(self):
        print "NO UNIT-TESTS FOR THIS MODULE"


######### Main function code #########

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    # taken:     --bBcCdDe---gGhH--------m-n-o---q-r-sStTuUvV--xX----  
    # free:      aB-------EfF----iIjJkKlL-M-N-OpP-Q-R--------wW--yYzZ  

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
    parser.add_option('-u', '--treat_unknown_as_match', action="store_true", default=False, 
                      help="When counting perfect reads, treat undefined alignment regions as matches (default %default)")
    parser.add_option('-U', '--dont_treat_unknown_as_match', action="store_false", dest='treat_unknown_as_match',
                      help="Turn -u off.")
    # LATER-TODO add a check to print a warning if any mutants are closer than X bases to each other - this seems to happen a lot with what looks like mutations, with the alignment locations just 1bp apart!  Should catch those and merge them somenow, or at least mark them as iffy in the output file?
    # MAYBE-TODO add a check to print a warning if any mutant has fewer than X% perfect reads; optionally mark/omit those mutants in the output file?
    # MAYBE-TODO add user-provided mutation cutoffs like in old deepseq_count_alignments.py, instead of just all reads and perfet reads?   parser.add_option('-m', '--mutation_cutoffs', default="1,3,10", metavar="<comma-separated-int-list>")

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

    ### gene-finding options 
    parser.add_option('-g', '--gene_position_reference_file', default=None, metavar='FILE', 
                      help="File to use to look up gene IDs based on chromosomal location (default %default)")
    parser.add_option('-G', '--gene_info_reference_file', default=None, metavar='FILE', 
                      help="File to use to look up gene names/descriptions from ID symbols (default %default)")
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

    ### output format options
    parser.add_option('-H', '--header_level', choices=['0','1','2'], default='2', metavar='0|1|2', 
                      help="Outfile header type:  0 - no header at all, 1 - a single line giving column headers, "
                          + "2 - full header with command, options, date, user etc (default %default) (also see -s)")
    parser.add_option('-n', '--N_sequences_per_group', type='int', default=2, metavar='N', 
                      help="How many most common sequences should be shown per group? (default %default)")
    parser.add_option('-s', '--add_summary_to_file', action="store_true", default=True, 
                      help="Print summary at the top of the file (default %default) (also see -H)")
    parser.add_option('-S', '--dont_add_summary_to_file', action="store_false", dest='add_summary_to_file', 
                      help="Turn -s off.")
    parser.add_option('-o', '--sort_data_key', choices=['position','read_count','none'], default='position', 
                      metavar='position|read_count|none', help="Sort the output data: by alignment position, read count, "
                             +"or don't sort at all (default %default) - sorting may be slow for large lists!")
    parser.add_option('-b', '--bad_chromosomes_count_only', default='insertion_cassette', metavar='comma-separated-list', 
                      help="Count reads aligning to these chromosomes and print the count in the header; "
                          +"otherwise treat them normally. (default %default) (also see -B)")
    parser.add_option('-B', '--bad_chromosomes_count_and_ignore', default='', metavar='comma-separated-list', 
                      help="Count reads aligning to these chromosomes and print the count in the header; "
                          +"otherwise ignore them and don't add to normal output. (default %default) (also see -b)")
    # LATER-TODO once I have a line-per-gene output format as well as a line-per-mutant one (separate dictionary/view in Insertional_mutant_library_dataset in mutant_analysis_classes.py, and separate output file):  Add extra options for that: count only mutants that are sense/antisense, only mutants in the exons/introns/UTRs, don't count mutants in the first/last X%/Xbp of the gene, do count mutants flanking the gene...

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
    
    

def run_main_function(infiles, outfile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """

    ### check/modify a few options
    # parse the -b/-B options
    bad_chromosomes_to_ignore = set(options.bad_chromosomes_count_and_ignore.split(',')) - set([''])
    bad_chromosomes_to_count = set(options.bad_chromosomes_count_only.split(',')) - set([''])

    ### generate empty alignment set object with basic read position/orientation properties defined by options
    all_alignment_data = mutant_analysis_classes.Insertional_mutant_library_dataset(options.read_cassette_end, 
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
        all_alignment_data.add_alignment_reader_to_data(infile_reader, options.input_collapsed_to_unique, 
                                     options.treat_unknown_as_match, bad_chromosomes_to_count, bad_chromosomes_to_ignore)

    ### optionally parse gene position/info files and look up the genes for each mutant in the data
    if options.gene_position_reference_file is not None:
        genefile = options.gene_position_reference_file
        if options.verbosity_level>1: print "adding genes from file %s to mutant data - time %s."%(genefile, time.ctime())
        all_alignment_data.find_genes_for_mutants(genefile, detailed_features=options.detailed_gene_features, 
                                                  known_bad_chromosomes=bad_chromosomes_to_count, 
                                                  N_run_groups=options.N_detail_run_groups, 
                                                  verbosity_level=options.verbosity_level)
        if options.gene_info_reference_file is not None:
            if options.verbosity_level>1: 
                print "adding gene details from file %s - time %s."%(options.gene_info_reference_file, time.ctime())
            all_alignment_data.add_detailed_gene_info(options.gene_info_reference_file)

    ### output
    # print summary info to stdout
    if options.verbosity_level>1:   print "\nDATA SUMMARY:"
    if options.verbosity_level>0:   all_alignment_data.print_summary()
    # print full data to outfile
    if options.verbosity_level>1:   print "printing output - time %s."%time.ctime()
    options.header_level = int(options.header_level)
    with open(outfile,'w') as OUTFILE:
        if options.header_level==2:
            write_header_data(OUTFILE,options)
        if options.add_summary_to_file:
            OUTFILE.write("### SUMMARY:\n")
            all_alignment_data.print_summary(OUTFILE, "#  ", "## ")
        if options.header_level==2:
            OUTFILE.write("### HEADER AND DATA:\n")
        elif options.add_summary_to_file:
            OUTFILE.write("### DATA:\n")
        header_line = True if options.header_level else False
        header_prefix = '' if options.header_level==1 else '# '
        all_alignment_data.print_data(OUTPUT=OUTFILE, sort_data_by=options.sort_data_key, 
                                      N_sequences=options.N_sequences_per_group, 
                                      header_line=header_line, header_prefix=header_prefix)


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

    run_main_function(infiles, outfile, options)

