#!/usr/bin/env python
""" Look at a gff file, print some overview info. 
 - print the basic structure of the records in the file (gene->mRNA->UTRs/exons or such)
 - print the gff IDs, sources, types etc and their counts
 - print the lengths (approximate) of each sequence, the number of genes it contains,
    and what percentage of the length is covered by genes
 - ...

WORK IN PROGRESS.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: gff_check_issues_coverage.py [options] gff_infile """

import sys, time
import unittest
import itertools
import pprint
from BCBio import GFF


def check_gene_coverage(sequence_record, length=None, check_for_overlap=True):
    """ Given a GFF sequence record (from BCBio.GFF parser), determine what fraction of its length is covered by genes.
    Length can be passed separately, or determined from the record (Caution: approximate! First-last feature only).
    """
    if length is None:  length = sequence_record.seq._length

    gene_length_total = 0
    for gene in sequence_record.features:
        gene_length_total += gene.location.end.position - gene.location.start.position
    gene_coverage_fraction = float(gene_length_total)/length

    # Check for overlapping genes and print a warning, since overlapping genes will make the measurement inaccurate
    if check_for_overlap:
        all_gene_positions = []
        for gene in sequence_record.features:
            all_gene_positions.append((gene.location.start.position,gene.location.end.position,gene.id))
    all_gene_positions.sort()
    for gene1pos,gene2pos in itertools.izip(all_gene_positions,all_gene_positions[1:]):
        if gene1pos[1]>=gene2pos[0]:
            print "WARNING: there are overlapping genes! (%s,%s) - %% of length covered by genes may not be accurate."\
                    %(gene1pos[2],gene2pos[2])
    # MAYBE-TODO actually adjust the measurement for overlapping genes?

    return gene_coverage_fraction 


######### Test code #########

class Testing_(unittest.TestCase):
    """ Unit-tests for _____. """
    def test__(self):
        pass
    # TODO implement!


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
    # TODO Not actually implemented!!  Just a copy from old code for reference.


######### Main function code #########

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as usage."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference. "
                          +"Ignores all other options/arguments. (default %default).")
    parser.add_option('-r', '--record_structure', action="store_true", default=False) 
    parser.add_option('-R', '--no_record_structure', action="store_false", dest='record_structure')
    parser.add_option('-c', '--type_counts', action="store_true", default=True) 
    parser.add_option('-C', '--no_type_counts', action="store_false", dest='type_counts')
    parser.add_option('-s', '--sequence_counts', action="store_true", default=False) 
    parser.add_option('-S', '--no_sequence_counts', action="store_false", dest='sequence_counts')
    parser.add_option('-o', '--source_counts', action="store_true", default=False) 
    parser.add_option('-O', '--no_source_counts', action="store_false", dest='source_counts')
    parser.add_option('-l', '--all_gff_limits', action="store_true", default=False) 
    parser.add_option('-L', '--no_all_gff_limits', action="store_false", dest='all_gff_limits')
    parser.add_option('-a', '--all_sequences', action="store_true", default=True) 
    parser.add_option('-A', '--no_all_sequences', action="store_false", dest='all_sequences')
    parser.add_option('-d', '--print_seq_details', action="store_true", default=False) 
    parser.add_option('-D', '--no_print_seq_details', action="store_false", dest='print_seq_details')

    return parser


def run_main_function(infile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument is generated by an optparse parser.
    """
    examiner = GFF.GFFExaminer()
    # I'm not sure why I need to open the file separately for each operation, but it doesn't work otherwise...

    if options.record_structure:
        with open(infile) as INFILE:
            print "\n *** Record structures ***"
            pprint.pprint(examiner.parent_child_map(INFILE))

    if options.type_counts:
        with open(infile) as INFILE:
            print "\n *** Type counts ***"
            pprint.pprint({'gff_type': examiner.available_limits(INFILE)['gff_type']})

    if options.sequence_counts:
        with open(infile) as INFILE:
            print "\n *** Sequence counts ***"
            pprint.pprint({'gff_id': examiner.available_limits(INFILE)['gff_id']})

    if options.source_counts:
        print "\n *** Source and source/type counts ***"
        with open(infile) as INFILE:
            pprint.pprint({'gff_source': examiner.available_limits(INFILE)['gff_source']})
        with open(infile) as INFILE:
            pprint.pprint({'gff_source_type': examiner.available_limits(INFILE)['gff_source_type']})

    if options.all_gff_limits:
        with open(infile) as INFILE:
            print "\n *** All GFF file limit field values ***"
            pprint.pprint(examiner.available_limits(INFILE))

    if options.all_sequences:
        with open(infile) as INFILE:
            print "\n *** Data for all sequences - Caution: looking at genes only! ***"
            print "       Caution: approximate sequence length is calculated from first to last gene only!)"
            genefile_parsing_limits = {'gff_type': ['gene']}
            for record in GFF.parse(INFILE, limit_info=genefile_parsing_limits):
                gene_coverage_percent = "%.0f%%"%(100*check_gene_coverage(record,None,True))
                print " * sequence %s: %s genes,"%(record.id, len(record.features)), 
                print "approximate length %s bp, with %s covered by genes"%(record.seq._length, gene_coverage_percent)
                if options.print_seq_details:   print record
                # MAYBE-TODO get the length from an actual genome fasta file - I think you can add those with GFF parser?

    # TODO check for weird cases! genes with more UTRs than expected, with no UTRs, with CDS on the outside of UTRs, ...  See notes_on_GFF_parsing.txt for some examples of this.
    # TODO make sure the gene IDs are all unique! 
    # TODO check for overlapping genes
    # TODO calculate the lowest distance between two genes

    # MAYBE-TODO there are probably going to be various alternative splicing issues...


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # if ran with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise parse the arguments and run main function
    try:                
        [infile] = args
    except ValueError:  
        parser.print_help()
        sys.exit("\nError: exactly one input gff file required!")

    run_main_function(infile, options)


