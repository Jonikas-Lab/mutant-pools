#!/usr/bin/env python
""" Convert a SAM file to a file containing genomic position, sequence, and total read count. 

The input file is a SAM deepseq alignment file created by bowtie, novoalign, or
other deepseq aligner programs (tested mainly on bowtie).  Multiple SAM input
files can be provided.

The output file is a tab-separated file, with one line per unique genomic
alignment location of the 3' end of the sequence (by default; other position
options may be used).  Each line will contain the following fields: chromosome,
position, most common non-mutated sequence, total number of aligned reads,
number of perfectly aligned reads.  (More output fields or separate files with
details on mutations, lengths, etc, may be added later.)

The program assumes the SAM file contains only unique matches (i.e. each read
was reported as aligning to at most one genomic location).
It also only deals with single-end alignments at the moment.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: deepseq_count_alignments.py [options] infile [infile2 infile3 ...] outfile """

# MAYBE-TODO add an option to make output go to STDOUT?
# MAYBE-TODO add mutation statistics and user-provided cutoffs like in the original old_deepseq_count_alignments.py script?

import sys, time
import unittest
from collections import defaultdict
from general_utilities import keybased_defaultdict, write_header_data
import HTSeq

######### NOTES ON THE SAM FORMAT
### Header:
# MAYBE-TODO do something useful with the SAM header?  Or at least copy it to outfile?
### Alignment line fields:
# * query template name
# * bitwise flag (relevant bits: 4 = unmapped, 16 = reverse-complement, 512 = failed quality control)
# * reference sequence name (* = unmapped)
# * leftmost mapping position (1-based) (0 = unmapped)
# * mapping quality ("-10 * log10(probability that position is wrong)"; 255 = unknown)
# * CIGAR string - descriptions of alignment matches/mismatches/etc (M/= match, I/D ins/del, X mismatch, S/H clipping)
# * (PE only - reference name of the mate fragment)
# * (PE only - position of the mate fragment)
# * template length
# * fragment sequence
# * ASCII of Phred-scaled base quality + 33   (original deepseq read quality)
# * OPTIONAL FIELDS, lots of different possibilities:   MD is mismatch info string, NM is edit distance to reference
#       (for info on MD field format see SAM manual footnote, and 
#        sam_MD_field_examples_*.txt files in experiments/reference_data/aligner_format_info)
#########

CIGAR_TYPES_MATCH = ['=']
CIGAR_TYPES_NOOP = ['S','H','P']
CIGAR_TYPES_MUTATION = ['X','I','D']
CIGAR_TYPES_INTRON = ['N']     # 'N' is for introns, but we shouldn't be paying attention to those for genomic DNA seq
CIGAR_TYPES_UNKNOWN = ['M']

VALID_POSITION_TYPES = ['leftmost','rightmost','5prime','3prime']

def check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as='unknown', ignore_introns=False):
    """ Return number of mutations in HTSeq_alignment, based on CIGAR string; -1 if unknown ('M') by default.
    If treat_unknown_as is 'unknown', return -1 whenever an unknown (M, may be match or mismatch) operation is found; 
     if treat_unknown_as is 'mutation' or 'match', count unknowns accordingly.  Return -1 if read is unaligned.
    If ignore_introns is False, count introns (N) as mutations; otherwise don't."""
    global CIGAR_TYPES_MUTATION, CIGAR_TYPES_INTRON, CIGAR_TYPES_UNKNOWN
    # just return -1 for unaligned reads
    if HTSeq_alignment.cigar is None:
        return -1
    # figure out whether to consider intron-skipping ('N') as a mutation or not, based on argument
    if ignore_introns:
        cigar_types_mutation = CIGAR_TYPES_MUTATION+CIGAR_TYPES_INTRON
    else:
        cigar_types_mutation = CIGAR_TYPES_MUTATION
    # figure out how to treat unknown matches ('M'), based on argument
    if treat_unknown_as=='unknown':
        cigar_types_unknown = CIGAR_TYPES_UNKNOWN
    elif treat_unknown_as=='mutation':
        cigar_types_mutation += CIGAR_TYPES_UNKNOWN
        cigar_types_unknown = []
    elif treat_unknown_as=='match':
        cigar_types_unknown = []
    else:
        raise ValueError("treat_unknown_as argument value must be 'mutation', 'match' or 'unknown'")
    # count the mutations, return total count (or instantly return -1 on finding an unknonw)
    mutations = 0
    for cigar_op in HTSeq_alignment.cigar:
        if cigar_op.type in cigar_types_mutation:
            mutations += cigar_op.size
        # if there's an unknown, just return -1, no need to count
        elif cigar_op.type in cigar_types_unknown:
            return -1
    return mutations


def check_mutation_count_by_optional_NM_field(HTSeq_alignment):
    """ Return number of mutations in HTSeq_alignment, based on optional NM field; -1 if unknown (NM field missing)."""
    # for unalign reads NM field is missing - returns -1
    try:                return HTSeq_alignment.optional_field('NM')
    except KeyError:    return -1


def check_mutation_count_by_optional_MD_field(HTSeq_alignment):
    """ Return number of mutations in HTSeq_alignment, based on optional MD field; -1 if unknown (MD field missing)."""
    # for info on MD field format see SAM manual footnote, 
    #   and sam_MD_field_examples_*.txt files in experiments/reference_data/aligner_format_info
    try:                mutation_string = HTSeq_alignment.optional_field('MD')
    except KeyError:    return -1
    # for unalign reads MD field is missing - returns -1
    mutation_letters = [c for c in mutation_string if not (c.isdigit() or c=='^')]
    return len(mutation_letters)


def check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as='unknown', ignore_introns=False):
    """ Return number of mutations in HTSeq_alignment (look at CIGAR string and NM and MD optional fields); -1 if unknown.
    First check the CIGAR string but only accept the answer if there are no unknown ('M') characters; 
     then check the NM and MD fields and return the result if those fields exist.
    If the CIGAR string is ambiguous and neither of the optional fields exist:
     - if treat_unknown_as is 'unknown', return -1
     - if treat_unknown_as is 'mutation' or 'match', return the CIGAR string result with unknowns counted accordingly.
    If ignore_introns is False, count introns (N) in CIGAR string as mutations; otherwise don't."""
    mutation_count = check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as='unknown', 
                                                          ignore_introns=ignore_introns)
    if not mutation_count==-1:  
        return mutation_count
    mutation_count = check_mutation_count_by_optional_NM_field(HTSeq_alignment)
    if not mutation_count==-1:  
        return mutation_count
    mutation_count = check_mutation_count_by_optional_MD_field(HTSeq_alignment)
    if not mutation_count==-1:  
        return mutation_count
    if treat_unknown_as=='unknown':     
        return -1
    return check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as=treat_unknown_as, 
                                                          ignore_introns=ignore_introns)


def get_chrom_and_pos_from_HTSeq_position(HTSeq_pos, pos_type):
    """ Return a (chrom, pos) tuple based on a HTSeq.GenomicPosition object, with pos being the location of the leftmost, rightmost, 5prime, or 3prime end of the read, depending on the value of pos_type."""
    if HTSeq_pos is None:
        sys.exit("Invalid position %s passed! Need an HTSeq iv object. (If empty, maybe read wasn't aligned?)"%HTSeq_pos)
    chrom = HTSeq_pos.chrom
    if pos_type=='leftmost':       pos = HTSeq_pos.start
    elif pos_type=='rightmost':       pos = HTSeq_pos.end
    elif pos_type=='5prime':    pos = HTSeq_pos.start if HTSeq_pos.strand=='+' else HTSeq_pos.end
    elif pos_type=='3prime':    pos = HTSeq_pos.end if HTSeq_pos.strand=='+' else HTSeq_pos.start
    else:                       raise ValueError("pos_type argument must be one of %s."%VALID_POSITION_TYPES)
    return chrom, pos


class Alignment_position_sequence_group():
    """ Data regarding sequences aligned to a particular genomic position (genomic position is set at initialization). 
    Variables: chromosome, position, total_read_count, perfect_read_count, sequence_counts (a sequence:count dictionary).
    Methods: add_read to add a given HTSeq read to the counts (doesn't check chromosome/position), 
     get_main_sequence to get the most common sequence from sequence_counts. """

    def __init__(self, chromosome, position):
        """ Set chromosome and position; initialize total_read_count, perfect_read_count and sequence_counts to 0 or {}."""
        self.chromosome, self.position = chromosome, position
        self.total_read_count = 0
        self.perfect_read_count = 0
        self.sequence_counts = defaultdict(lambda: 0)

    def add_read(self, HTSeq_alignment, treat_unknown_as_match=False):
        """ Add a read to the data: increment total_read_count, increment perfect_read_count if read is a perfect 
        alignment, increment the appropriate field of sequence_counts based on read sequence."""
        self.sequence_counts[HTSeq_alignment.read.seq] += 1
        self.total_read_count += 1
        # figure out whether the read is perfect, treating unknowns ('M' in CIGAR string) as desired
        treat_unknown_as = 'match' if treat_unknown_as_match else 'mutation'
        mutation_count = check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as=treat_unknown_as)
        if mutation_count==0:  
            self.perfect_read_count += 1

    def get_main_sequence(self, N=1):
        """ Return the most common sequence in this group and its count (or Nth most common sequence if N is provided)."""
        sequences_by_count = sorted([(count,seq) for (seq,count) in self.sequence_counts.iteritems()], reverse=True)
        # try returning the Nth sequence and count; return nothing if there are under N sequences.
        try:                return (sequences_by_count[N-1][1], sequences_by_count[N-1][0])
        except IndexError:  return ('',0)


class All_alignments_grouped_by_pos():
    """ Essentially a dictionary of alignment_position_sequence_group with position data (chrom,pos) as keys. """

    def __init__(self, position_type):
        """ Checks position_type and assigns to self.position_type; initializes self.alignment_position_data_dict.
        position_type must be one of %s.
        self.alignment_position_data_dict is a dictionary that generates a new alignment_position_sequence_group object
         based on the key (i.e. position) if the key isn't already in the dictionary.  """%VALID_POSITION_TYPES
        self.alignment_position_data_dict = keybased_defaultdict(lambda key: Alignment_position_sequence_group(*key))
        if not position_type in VALID_POSITION_TYPES: 
            raise ValueError("The position_type variable must be one of %s!"%VALID_POSITION_TYPES)
        self.position_type = position_type
        self.total_read_count, self.aligned_read_count, self.unaligned_read_count = 0,0,0

    def add_alignment_reader_to_data(self, HTSeq_alignment_reader, treat_unknown_as_match=False):
        """ Adds all alignments to self.alignment_position_data_dict based on position.
        Input must be a list/generator/etc of HTSeq.Alignment objects (usually an HTSeq.SAM_Reader)."""
        for aln in HTSeq_alignment_reader:
            self.total_read_count += 1
            if (not aln.aligned) or (aln.iv is None):
                self.unaligned_read_count += 1
                continue
            self.aligned_read_count += 1
            position_info = get_chrom_and_pos_from_HTSeq_position(aln.iv, self.position_type)
            self.alignment_position_data_dict[position_info].add_read(aln, treat_unknown_as_match=treat_unknown_as_match)

    def print_summary(self, OUTPUT=None, line_prefix = ''):
        """ Print basic read and group counts (prints to stdout by default, can also pass an open file object)."""
        if OUTPUT is None:
            OUTPUT = sys.stdout
        OUTPUT.write("%s Total reads processed: %s\n"%(line_prefix, self.total_read_count))
        OUTPUT.write("%s Aligned reads: %s\n"%(line_prefix, self.aligned_read_count))
        OUTPUT.write("%s Unaligned reads: %s\n"%(line_prefix, self.unaligned_read_count))
        position_info = "looking at %s of read"%self.position_type
        OUTPUT.write("%s Read groups by alignment position (%s): %s\n"%(line_prefix, position_info, 
                                                                        len(self.alignment_position_data_dict)))

    def print_data(self, OUTPUT=None, N_sequences=2, header_line=True, header_prefix="# "):
        """ Print the full data:  the read count for each group of sequences with the same position.
        If N_sequences<1, only prints position (chromosome and position), total and perfect read count, 
          and the number of sequence variants.  If N_sequences==1, also prints the most common sequence and count; 
          if N_sequences>1, adds the second most common sequence and count, and so on.
        Output is tab-separated, with headers starting with "# ".  Prints to an open file object (stdout by default).
        """
        if OUTPUT is None:
            OUTPUT = sys.stdout

        if header_line:
            headers = ['chromosome','position','total_reads','perfect_reads', 'N_sequence_variants']
            for N in range(1,N_sequences+1):
                headers.extend(['sequence_%s'%N,'count_seq_%s'%N])
            OUTPUT.write(header_prefix + '\t'.join(headers) + "\n")

        for group in self.alignment_position_data_dict.itervalues():
            group_data = [group.chromosome, group.position, group.total_read_count, 
                          group.perfect_read_count, len(group.sequence_counts)]
            OUTPUT.write('\t'.join([str(x) for x in group_data]))
            for N in range(1,N_sequences+1):
                OUTPUT.write('\t%s\t%s'%group.get_main_sequence(N))
            OUTPUT.write('\n')


def run_main_function(infiles, outfile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument is generated by an optparse parser.
    """
    all_alignment_data = All_alignments_grouped_by_pos(options.position_type)

    for infile in infiles:
        if options.verbose: print "parsing input file %s - time %s."%(infile, time.ctime())
        infile_reader = HTSeq.SAM_Reader(infile)
        all_alignment_data.add_alignment_reader_to_data(infile_reader, options.treat_unknown_as_match)
        if options.verbose: print "finished parsing input file %s - time %s."%(infile, time.ctime())

    all_alignment_data.print_summary()

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
        all_alignment_data.print_data(OUTFILE, options.N_sequences_per_group, header_line, header_prefix)


### Test code

class Testing_single_functions(unittest.TestCase):
    """ Unit-tests for all the top-level functions in this module. """

    def test__check_mutation_count_by_CIGAR_string(self):
        pass
        # TODO implement!

    def test__check_mutation_count_by_optional_NM_field(self):
        pass
        # TODO implement!

    def test__check_mutation_count_by_optional_MD_field(self):
        pass
        # TODO implement!

    def test__check_mutation_count_try_all_methods(self):
        pass
        # TODO implement!

    def test__get_chrom_and_pos_from_HTSeq_position(self):
        pass
        # TODO implement!


class Testing_Alignment_position_sequence_group(unittest.TestCase):
    """ Unit-tests for the Alignment_position_sequence_group class and its methods. """

    pass
    # TODO implement!


class Testing_All_alignments_grouped_by_pos(unittest.TestCase):
    """ Unit-tests for the All_alignments_grouped_by_pos class and its methods. """

    pass
    # TODO implement!


# TODO Write overall test comparing to current output or to known small outputs?

if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    # optparse only used here!  The options are converted to normal arguments before passing to main.
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default False).")
    parser.add_option('-p', '--position_type', choices=VALID_POSITION_TYPES, 
                      default='3prime', metavar='|'.join(VALID_POSITION_TYPES), 
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
    parser.add_option('-u', '--treat_unknown_as_match', action="store_true", default=False, 
                      help="When counting perfect reads, treat undefined alignment regions as matches (default %default)")
    parser.add_option('-U', '--dont_treat_unknown_as_match', action="store_false", dest='treat_unknown_as_match',
                      help="Turn -u off.")
    # TODO add some way of specifying chromosomes or chromosome regions to ignore?  Like insertion_cassette
    # TODO separator/format in case I want csv files instead of tab-separated ones?  I'd like to be able to read this file by eye, so I'd like to be able to have a comma-separated format with arbitrary whitespace to line things up.
    parser.add_option('-v', '--verbose', action="store_true", default=False)
    (options, args) = parser.parse_args()

    # if ran with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in simple test suite.")
        # tun unittest.main, passing it no arguments (by default it takes sys.argv and complains about options)
        unittest.main(argv=[sys.argv[0]])


    # otherwise parse the arguments and run main function
    if len(args)<2:    
        parser.print_help()
        sys.exit("\nError: at least one infile and exactly one outfile are required!")
    infiles = args[:-1]
    outfile = args[-1]

    run_main_function(infiles, outfile, options)


    # TODO we DON'T want the grouping to be based only on alignment location, I think!  What if we have two mutants that inserted into the same location but in opposite orientations?  12--->345 and 123<---45 - position should be "3" in both cases (is that really how it comes out? check!), but I think those are separate mutants and should be treated separately, right?  Unlikely, of course, but still.

    ### MAYBE-TODO more options I might want (see old_deepseq_count_alignments.py for old code for dealing with them)
    #parser.add_option('-m', '--mutation_cutoffs', default="1,3,10", metavar="<comma-separated-int-list>")
    #parser.add_option('-r', '--count_repeat_matches', choices=["none","all","by_seq","from_file"], default="none", metavar="[none|all|by_seq|from_file]", help="How to deal with cases of one read aligning to multiple reference sequences: 'all' - count it as matching both sequences, 'none' - ignore it, 'by_seq' - only count it if both reference sequences are identical, 'from_file' - only count it if the pair occurs in the allowed_repeats_file (see -a option).")
    #parser.add_option('-a', '--allowed_repeats_file', type="string", metavar="FILE", help="See -r option 'from_file' for the explanation. The file should contain any number of lines of the format \"sequence<tab> name1,name2,name3,...\", and may be generated by the find_library_duplicates.py script.")
    # MAYBE-TODO alignment filetype option? (SAM or bowtie or novoalign - currently only implementing SAM)
