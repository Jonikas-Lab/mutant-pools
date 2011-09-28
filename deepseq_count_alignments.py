#!/usr/bin/env python
""" Convert a SAM file to a file containing genomic position, total read count,
and optionally other details. 

The input file is a SAM deepseq alignment file created by bowtie, novoalign, or
other deepseq aligner programs (tested mainly on bowtie).  Multiple SAM input
files can be provided.

The output file is a tab-separated file, with one line per unique genomic
alignment location of the 3' end of the sequence (by default; other position
options may be used).  Each line will contain the following fields: chromosome,
position, total number of aligned reads, number of perfectly aligned reads.
Optionally: most common sequence and its count, second most common sequence and
its count, etc.  (More output fields or separate files with details on
mutations, lengths, etc, may be added later.)

The program assumes the SAM file contains only unique matches (i.e. each read
was reported as aligning to at most one genomic location).
It also only deals with single-end alignments at the moment.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: deepseq_count_alignments.py [options] infile [infile2 infile3 ...] outfile """

import sys, os, time
import filecmp
import unittest
from general_utilities import write_header_data
import HTSeq
import deepseq_analysis_classes


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    test_runs = [("-H 1 -s -n3 -o -q -U -p leftmost","test_data/test_input.sam","test_data/test_output__U_leftmost.txt"),
                 ("-H 1 -s -n3 -o -q -U -p rightmost","test_data/test_input.sam","test_data/test_output__U_rightmost.txt"),
                 ("-H 1 -s -n3 -o -q -U -p 5prime","test_data/test_input.sam","test_data/test_output__U_5prime.txt"),
                 ("-H 1 -s -n3 -o -q -U -p 3prime","test_data/test_input.sam","test_data/test_output__U_3prime.txt"),
                 ("-H 1 -s -n3 -o -q -u -p leftmost","test_data/test_input.sam","test_data/test_output__u_leftmost.txt")]
    #  (using -s and -H 1 to get all relevant info but not have to deal with changing timestamps/etc)
    for option_string, infile, reference_file in test_runs:
        print(" * New test run, with options: %s (infile %s, reference outfile %s)"%(option_string,infile,reference_file))
        # regenerate options with test argument string
        parser = define_option_parser()
        (options, _) = parser.parse_args(option_string.split())
        outfile = "test_data/test_output.txt"
        run_main_function([infile], outfile, options)
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
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                      + "Ignores all other options/arguments. (default %default).")
    parser.add_option('-p', '--position_type', choices=deepseq_analysis_classes.VALID_POSITION_TYPES, default='3prime', 
                      metavar='|'.join(deepseq_analysis_classes.VALID_POSITION_TYPES), 
                      help="Which position feature should be used to group reads together? (default %default) "
                           + "leftmost/rightmost refer to where the first aligned base of the read lies on the reference, "
                           + "regardless of read orientation; 5prime/3prime is by position of specific end of the read.")
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
    parser.add_option('-u', '--treat_unknown_as_match', action="store_true", default=False, 
                      help="When counting perfect reads, treat undefined alignment regions as matches (default %default)")
    parser.add_option('-U', '--dont_treat_unknown_as_match', action="store_false", dest='treat_unknown_as_match',
                      help="Turn -u off.")
    # TODO if the imput data had been ran through fastx_collapser before alignment, each unique sequence will have exactly one count - need extract the actual original read count data from the read name!!  Implement in deepseq_analysis_classes.py Alignment_position_sequence_group.add_read
    # TODO add some way of specifying chromosomes or chromosome regions to ignore?  Like insertion_cassette
    # TODO separator/format in case I want csv files instead of tab-separated ones?  I'd like to be able to read this file by eye, so I'd like to be able to have a comma-separated format with arbitrary whitespace to line things up.
    # MAYBE-TODO add user-provided mutation cutoffs like in old_deepseq_count_alignments.py, instead of just all reads and perfet reads
    # MAYBE-TODO add a check to print a warning if any mutants are closer than X bases to each other; optionally mark/omit those mutants in the output file?
    # MAYBE-TODO add a check to print a warning if any mutant has fewer than X% perfect reads; optionally mark/omit those mutants in the output file?
    # TODO eventually I want to implement grouping based on sequence (clustering?) instead of just based on alignment position!  See "Notes on grouping mutants based on sequence/position/etc" section in ../notes.txt
    parser.add_option('-q', '--quiet', action="store_true", default=False, help="Don't print summary to STDOUT.")
    parser.add_option('-v', '--verbose', action="store_true", default=False, help="Print progress reports to STDOUT.")

    return parser


def run_main_function(infiles, outfile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument is generated by an optparse parser.
    """
    all_alignment_data = deepseq_analysis_classes.All_alignments_grouped_by_pos(options.position_type)

    for infile in infiles:
        if options.verbose: print "parsing input file %s - time %s."%(infile, time.ctime())
        infile_reader = HTSeq.SAM_Reader(infile)
        all_alignment_data.add_alignment_reader_to_data(infile_reader, options.treat_unknown_as_match)
        if options.verbose: print "finished parsing input file %s - time %s."%(infile, time.ctime())

    if not options.quiet:   all_alignment_data.print_summary()

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
    if len(args)<2:    
        parser.print_help()
        sys.exit("\nError: at least one infile and exactly one outfile are required!")
    infiles = args[:-1]
    outfile = args[-1]

    run_main_function(infiles, outfile, options)


    # TODO we DON'T want the grouping to be based only on alignment location, I think!  What if we have two mutants that inserted into the same location but in opposite orientations?  12--->345 and 123<---45 - position should be "3" in both cases (is that really how it comes out? check!), but I think those are separate mutants and should be treated separately, right?  Unlikely, of course, but still.

    # MAYBE-TODO do I want to try outputting this in some existing bioinformatics format instead of a made-up one?

    # MAYBE-TODO add an option to make output go to STDOUT?

    ### MAYBE-TODO more options I might want (see old_deepseq_count_alignments.py for old code for dealing with them)
    #parser.add_option('-m', '--mutation_cutoffs', default="1,3,10", metavar="<comma-separated-int-list>")
    #parser.add_option('-r', '--count_repeat_matches', choices=["none","all","by_seq","from_file"], default="none", metavar="[none|all|by_seq|from_file]", help="How to deal with cases of one read aligning to multiple reference sequences: 'all' - count it as matching both sequences, 'none' - ignore it, 'by_seq' - only count it if both reference sequences are identical, 'from_file' - only count it if the pair occurs in the allowed_repeats_file (see -a option).")
    #parser.add_option('-a', '--allowed_repeats_file', type="string", metavar="FILE", help="See -r option 'from_file' for the explanation. The file should contain any number of lines of the format \"sequence<tab> name1,name2,name3,...\", and may be generated by the find_library_duplicates.py script.")
    # MAYBE-TODO alignment filetype option? (SAM or bowtie or novoalign - currently only implementing SAM)
