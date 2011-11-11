#!/usr/bin/env python
""" Group the reads in a deepseq alignment file by genomic alignment position;
output a new file with one line per group, with total read counts, most common
sequence variants, and optionally other details. 

The input file is a SAM deepseq alignment file created by bowtie, novoalign, or
other deepseq aligner programs (tested mainly on bowtie).

The output file is a tab-separated file, with one line per unique genomic
alignment location of the 3' end of the sequence (by default; other position
options may be used).  Each line will contain the following fields: chromosome,
strand, position, gene, orientation vs gene, gene feature, total number of
aligned reads, number of perfectly aligned reads.  Optionally: most common
sequence and its count, second most common sequence and its count, etc.  (More
output fields or separate files with details on mutations, lengths, etc, may be
added later.) The gene information is only filled in if a gene annotation file
is provided via options.

The program assumes the SAM file contains only unique matches (i.e. each read
was reported as aligning to at most one genomic location).
It also only deals with single-end alignments at the moment.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: deepseq_count_alignments.py [options] infile outfile """

# basic libraries
import sys, os, time
import filecmp
import unittest
# other libraries
import HTSeq
# my modules
from general_utilities import write_header_data
import deepseq_analysis_classes


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    test_infile = "test_data/test_input.sam"
    test_constant_options = "-H 1 -s -n3 -o -q"
    #  (using -s and -H 1 to get all relevant info but not have to deal with changing timestamps/etc)
    test_runs = [("-U -p leftmost", "test_data/test_output__leftmost.txt"),
                 ("-U -p rightmost", "test_data/test_output__rightmost.txt"),
                 ("-U -p 5prime", "test_data/test_output__5prime.txt"),
                 ("-U -p 3prime", "test_data/test_output__3prime.txt"),
                 ("-u -p leftmost", "test_data/test_output__u.txt"),
                 ("-b special_chromosome -U -p leftmost", "test_data/test_output__b.txt"),
                 ("-B special_chromosome -U -p leftmost", "test_data/test_output__B.txt")]
    for variable_option_string, reference_file in test_runs:
        option_string = test_constant_options + ' ' + variable_option_string
        print(" * New test run, with options: %s (infile %s, reference outfile %s)"%(option_string,test_infile,reference_file))
        # regenerate options with test argument string
        parser = define_option_parser()
        (options, _) = parser.parse_args(option_string.split())
        outfile = "test_data/test_output.txt"
        run_main_function(test_infile, outfile, options)
        # compare outfile to reference file: remove outfile and keep going if correct, otherwise exit with message.
        if filecmp.cmp(outfile, reference_file, shallow=False):
            os.remove(outfile)
        else:
            print("TEST FAILED!!  Reference file %s and output file %s differ - PLEASE COMPARE."%(reference_file,outfile))
            return 1
    print("*** Test runs finished - EVERYTHING IS FINE. ***")
    return 0


######### Main function code #########

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as usage."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    ### functionality options
    parser.add_option('-p', '--position_type', choices=deepseq_analysis_classes.VALID_POSITION_TYPES, default='3prime', 
                      metavar='|'.join(deepseq_analysis_classes.VALID_POSITION_TYPES), 
                      help="Which position feature should be used to group reads together? (default %default) "
                          + "leftmost/rightmost refer to where the first aligned base of the read lies on the reference, "
                          + "regardless of read orientation; 5prime/3prime is by position of specific end of the read.")
    parser.add_option('-u', '--treat_unknown_as_match', action="store_true", default=False, 
                      help="When counting perfect reads, treat undefined alignment regions as matches (default %default)")
    parser.add_option('-U', '--dont_treat_unknown_as_match', action="store_false", dest='treat_unknown_as_match',
                      help="Turn -u off.")
    # TODO add a check to print a warning if any mutants are closer than X bases to each other - this seems to happen a lot with what looks like mutations, with the alignment locations just 1bp apart!  Should catch those and merge them somenow, or at least mark them as iffy in the output file?
    # MAYBE-TODO add a check to print a warning if any mutant has fewer than X% perfect reads; optionally mark/omit those mutants in the output file?
    # LATER-TODO eventually I want to implement grouping based on sequence (clustering?) instead of just based on alignment position!  See "Notes on grouping mutants based on sequence/position/etc" section in ../notes.txt

    ### input options
    parser.add_option('-c','--input_collapsed_to_unique', action='store_true', default=False, 
                      help="Use to get correct original total read counts if the data was collapsed to unique sequences using fastx_collapser before alignment (default %default).")
    parser.add_option('-C','--input_not_collapsed_to_unique', action='store_false', dest="input_collapsed_to_unique", 
                      help="Turn -c off.")
    parser.add_option('-m', '--input_metadata_file', default='AUTO', metavar='FILE', 
                      help="File containing preprocessing and alignment metadata (scripts/options used etc). "
                          +"Default: <infile_basename>_info.txt. Warning will be raised if not found. "
                          +"Pass NONE to not look for a metadata file at all.")

    ### gene-finding options 
    parser.add_option('-g', '--gene_position_reference_file', default=None, metavar='FILE', 
                      help="File to use to look up gene IDs based on chromosomal location (default %default)")
    parser.add_option('-G', '--gene_info_reference_file', default=None, metavar='FILE', 
                      help="File to use to look up gene names/descriptions from ID symbols (default %default)"
                          +"       NOT IMPLEMENTED")
    # LATER-TODO implement -G! Get gene name and hopefully description, and maybe GO factors and things...
    parser.add_option('-d', '--detailed_gene_features', action="store_true", default=False,
                      help="Find out what part of the gene (UTR,intron,exon) a mutant hit, based on the -g file "
                          +"- may take a LOT of memory! (default %default)"
                          +"      NOT IMPLEMENTED")
    # TODO implement -d! GFF parsing already works, the problem is just that it seems to take up a lot of memory... 
    #   also sometimes gene structure is CONFUSING, there are splice variants and weird-positioned UTRs/exons and what not.
    #   if I do get this, where should it go?  Possibly in parentheses next to gene name in standard by-mutant output, 
    #     and I suppose in by-gene output there should be a column for number of mutants hitting various features...
    parser.add_option('-D', '--no_detailed_gene_features', action="store_false", dest='detailed_gene_features',
                      help="Turns -d off.")
    # MAYBE-TODO add a "flank" option (with variable size), to catch mutants that are in the flanks of genes? Do we care?
    # MAYBE-TODO add a "negative flank" option (with variable size), to ignore mutants that are in the start/end of genes?

    ### output format options
    parser.add_option('-H', '--header_level', choices=['0','1','2'], default='2', metavar='0|1|2', 
                      help="Outfile header type:  0 - no header at all, 1 - a single line giving column headers, "
                          + "3 - full header with command, options, date, user etc (default %default) (also see -s)")
    parser.add_option('-n', '--N_sequences_per_group', type='int', default=2, metavar='N', 
                      help="How many most common sequences should be shown per group? (default %default)")
    parser.add_option('-s', '--add_summary_to_file', action="store_true", default=True, 
                      help="Print summary at the end of the file (default %default) (also see -H)")
    parser.add_option('-S', '--dont_add_summary_to_file', action="store_false", dest='add_summary_to_file', 
                      help="Turn -s off.")
    parser.add_option('-o', '--sort_data_by_position', action="store_true", default=False, 
                      help="Sort the output data by alignment position (default %default) - CAUTION: MAY BE SLOW!")
    parser.add_option('-O', '--dont_sort_data_by_position', action="store_false", dest='sort_data_by_position', 
                      help="Turn -o off.")
    parser.add_option('-b', '--bad_chromosomes_count_only', default='insertion_cassette', metavar='comma-separated-list', 
                      help="Count reads aligning to these chromosomes and print the count in the header; "
                          +"otherwise treat them normally. (default %default) (also see -B)")
    parser.add_option('-B', '--bad_chromosomes_count_and_ignore', default='', metavar='comma-separated-list', 
                      help="Count reads aligning to these chromosomes and print the count in the header; "
                          +"otherwise ignore them and don't add to normal output. (default %default) (also see -b)")
    # LATER-TODO should have a line-per-gene output format as well as a line-per-mutant one!  Which should be a separate dictionary/view in All_alignments_grouped_by_pos in deepseq_analysis_classes.py, I suppose.  Extra options for that: count only mutants that are sense/antisense, only mutants in the exons/introns/UTRs, don't count mutants in the first/last X%/Xbp of the gene, do count mutants flanking the gene...
        

    # MAYBE-TODO add user-provided mutation cutoffs like in old_deepseq_count_alignments.py, instead of just all reads and perfet reads?   parser.add_option('-m', '--mutation_cutoffs', default="1,3,10", metavar="<comma-separated-int-list>")
    parser.add_option('-q', '--quiet', action="store_true", default=False, help="Don't print summary to STDOUT.")
    parser.add_option('-v', '--verbose', action="store_true", default=False, help="Print progress reports to STDOUT.")

    return parser


def run_main_function(infile, outfile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """
    ### generate empty alignment set object
    all_alignment_data = deepseq_analysis_classes.All_alignments_grouped_by_pos(options.position_type)

    ### parse preprocessing/alignment metadata file to get discarded read count, pass it to all_alignment_data
    # all_alignment_data initializes it to 'unkown', so if file is not given or can't be found, no need to do anything
    if options.input_metadata_file is not 'NONE':
        if options.input_metadata_file == 'AUTO':
            options.input_metadata_file = os.path.splitext(infile)[0] + '_info.txt'
            if options.verbose:  print 'Automatically determining metadata input file name: %s'%options.input_metadata_file
        else:
            if options.verbose:  print 'Metadata input file name provided in options: %s'%options.input_metadata_file
        if not os.path.exists(options.input_metadata_file):
            print 'Warning: metadata input file %s not found! Proceeding without it.'%options.input_metadata_file
        else:
            for line in open(options.input_metadata_file):
                if line.startswith('## reads removed: '):
                    all_alignment_data.add_discarded_reads(line.split()[3])
                    break
            if all_alignment_data.discarded_read_count=='unknown':
                print 'Warning: discarded read count not found in metadata input file %s!'%options.input_metadata_file
    # MAYBE-TODO get the final total number of reads from the metadata infile and make sure it's the same 
    #   as the number of processed reads I get from all_alignment_data.print_summary()?
    
    ### parse input file and store data - the add_alignment_reader_to_data function here does pretty much all the work!
    # parse the -b/-B options
    bad_chromosomes_to_ignore = set(options.bad_chromosomes_count_and_ignore.split(',')) - set([''])
    bad_chromosomes_to_count = set(options.bad_chromosomes_count_only.split(',')) - set([''])
    # initialize a parser for the SAM infile
    if options.verbose: print "parsing input file %s - time %s."%(infile, time.ctime())
    infile_reader = HTSeq.SAM_Reader(infile)
    # fill the new alignment set object with data from the infile parser
    all_alignment_data.add_alignment_reader_to_data(infile_reader, options.input_collapsed_to_unique, 
                                     options.treat_unknown_as_match, bad_chromosomes_to_count, bad_chromosomes_to_ignore)
    
    ### optionally parse gene position/info files and look up the genes for each mutant in the data
    if options.gene_position_reference_file:
        genefile = options.gene_position_reference_file
        if options.verbose: print "parsing gene file %s and adding data from it - time %s."%(genefile, time.ctime())
        all_alignment_data.add_gene_positions_to_data(genefile, detailed_features=options.detailed_gene_features, 
                                                      gene_info_file=options.gene_info_reference_file, 
                                                      known_bad_chromosomes=bad_chromosomes_to_count)

    ### output
    # print summary info to stdout
    if not options.quiet:   
        all_alignment_data.print_summary()
    # print full data to outfile
    if options.verbose: print "printing output - time %s."%time.ctime()
    options.header_level = int(options.header_level)
    with open(outfile,'w') as OUTFILE:
        if options.header_level==2:
            write_header_data(OUTFILE,options)
        if options.add_summary_to_file:
            OUTFILE.write("### SUMMARY:\n")
            all_alignment_data.print_summary(OUTFILE, "# ")
        if options.header_level==2:
            OUTFILE.write("### HEADER AND DATA:\n")
        elif options.add_summary_to_file:
            OUTFILE.write("### DATA:\n")
        header_line = True if options.header_level else False
        header_prefix = '' if options.header_level==1 else '# '
        all_alignment_data.print_data(OUTPUT=OUTFILE, sort_data=options.sort_data_by_position, 
                                      N_sequences=options.N_sequences_per_group, 
                                      header_line=header_line, header_prefix=header_prefix)


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # if ran with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        print("      (testing both the deepseq_analysis_classes.py module and this module)")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(deepseq_analysis_classes)
        unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise parse the arguments and run main function
    try:
        [infile,outfile] = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: exactly one infile and exactly one outfile are required!")
    # MAYBE-TODO allow it to take multiple infiles if metafile is not provided?  Just in case I want to merge different samples for some reason or another.  Could do it if metafile is AUTO, too, really - I don't include metafile contents in output, just the one discarded read number, and I can add that from all the metafiles together.

    run_main_function(infile, outfile, options)

    # MAYBE-TODO If we have two mutants that inserted into the exact same location in opposite directions, WOULD they actually get reported as the same position, or offset by one?  Offset by one, probably: 12-->345 would deepseq "345" and get 3 as the position; 12345 would deepseq "12" and get 2 as the position, I think.  May want to fix that. Again, though, this is EXTREMELY UNLIKELY to make a difference (would have to be right on the edge of a gene/feature/something).
