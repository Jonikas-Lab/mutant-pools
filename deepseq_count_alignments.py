#!/usr/bin/env python
""" Convert a SAM file to a file containing genomic position, sequence, and total read count. 

The input file is a SAM deepseq alignment file created by bowtie, novoalign, or
other deepseq aligner programs (tested mainly on bowtie).  Multiple SAM input
files can be provided.

The output file is a tab-separated file, with one line per unique genomic
alignment location of the 3' end of the sequence (?).  Each line will contain
the following fields: chromosome, position, most common non-mutated sequence,
total number of aligned reads, number of perfectly aligned reads.  (More output
fields or separate files with details on mutations, lengths, etc, may be added
later.)

The program assumes the SAM file contains only unique matches (i.e. each read
was reported as aligning to at most one genomic location).
It also only deals with single-end alignments at the moment.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: deepseq_count_alignments.py [options] infile [infile2 infile3 ...] outfile """

# MAYBE-TODO add an option to make output go to STDOUT?
# MAYBE-TODO add mutation statistics and user-provided cutoffs like in the original old_deepseq_count_alignments.py script?

import HTSeq
from collections import defaultdict
from general_utilites import keybased_defaultdict

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

def check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as='unknown', ignore_introns=False):
    """ Return number of mutations in HTSeq_alignment, based on CIGAR string; -1 if unknown ('M') by default.
    If treat_unknown_as is 'unknown', return -1 whenever an unknown (M, may be match or mismatch) operation is found; 
     if treat_unknown_as is 'mutation' or 'match', count unknowns accordingly.
    If ignore_introns is False, count introns (N) as mutations; otherwise don't."""
    global CIGAR_TYPES_MUTATION, CIGAR_TYPES_INTRON, CIGAR_TYPES_UNKNOWN
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
        if cigar.op.type in cigar_types_mutation:
            mutations += cigar.op.size
        # if there's an unknown, just return -1, no need to count
        elif cigar.op.type in cigar_types_unknown:
            return -1
    return mutations
# TODO what does the CIGAR string of an unaligned read look like? 


def check_mutation_count_by_optional_NM_field(HTSeq_alignment):
    """ Return number of mutations in HTSeq_alignment, based on optional NM field; -1 if unknown (NM field missing)."""
    try:                return HTSeq_alignment.optional_field('NM')
    except KeyError:    return -1
# TODO what does the NM field of an unaligned read look like? 


def check_mutation_count_by_optional_MD_field(HTSeq_alignment):
    """ Return number of mutations in HTSeq_alignment, based on optional MD field; -1 if unknown (MD field missing)."""
    # for info on MD field format see SAM manual footnote, 
    #   and sam_MD_field_examples_*.txt files in experiments/reference_data/aligner_format_info
    try:                mutation_string = HTSeq_alignment.optional_field('MD')
    except KeyError:    return -1
    mutation_letters = [c for c in mutation_string if not (c.isdigit() or c=='^')]
    return len(mutation_letters)
# TODO what does the MD field of an unaligned read look like? 


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
    """ Return a (chrom, pos) tuple from an HTSeq.GenomicPosition object; pos_type can be start, end, 3prime, 5prime."""
    chrom = HTSeq_pos.chrom
    if pos_type=='start':       pos = HTSeq_pos.start
    elif pos_type=='end':       pos = HTSeq_pos.end
    elif pos_type=='5prime':    pos = HTSeq_pos.start if HTSeq_pos.strand=='+' else HTSeq_pos.end
    elif pos_type=='3prime':    pos = HTSeq_pos.end if HTSeq_pos.strand=='+' else HTSeq_pos.start
    else:                       raise ValueError("pos_type argument must be 'start', 'end', '3prime', or '5prime'.")
    return chrom, pos
# TODO what does the position of an unaligned read look like? 


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
        self.sequence_counts(HTSeq_alignment.read.seq) += 1
        self.total_read_count += 1
        # figure out whether the read is perfect, treating unknowns ('M' in CIGAR string) as desired
        treat_unknown_as = 'match' if treat_unknown_as_match else 'mutation'
        mutation_count = check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as=treat_unknown_as)
        if mutation_count==0:  
            self.perfect_read_count += 1

    def get_main_sequence(self):
        """ Return the most common sequence in this group."""
        sequences_by_count = sorted([(count,seq) for (seq,count) in self.sequence_counts.iteritems()], reverse=True)
        return sequences_by_count[0][1]


class All_alignments_grouped_by_pos():
    """ Essentially a dictionary of alignment_position_sequence_group with position data (chrom,pos) as keys. """

    def __init__(self, position_type):
        """ Checks position_type and assigns to self.position_type; initializes self.alignment_position_data_dict.
        position_type must be 'start','end','3prime',or '5prime'.
        self.alignment_position_data_dict is a dictionary that generates a new alignment_position_sequence_group object
         based on the key (i.e. position) if the key isn't already in the dictionary.  """
        self.alignment_position_data_dict = keybased_defaultdict(lambda key: alignment_position_sequence_group(*key))
        if not position_type in ['start','end','3prime','5prime']: 
            raise ValueError("The position_type variable must be 'start','end','3prime',or '5prime'!")
        self.position_type = position_type
        self.unaligned_count = 0

    def add_alignment_reader_to_data(self, HTSeq_alignment_reader, treat_unknown_as_match=False):
        """ Adds all alignments to self.alignment_position_data_dict based on position.
        Input must be a list/generator/etc of HTSeq.Alignment objects (usually an HTSeq.SAM_Reader)."""
        for aln in HTSeq_alignment_reader:
            if not aln.aligned:
                self.unaligned_count += 1
                continue
            position_info = get_chrom_and_pos_from_HTSeq_position(aln.iv, self.position_type)
            self.alignment_position_data_dict[position_info].add_read(aln, treat_unknown_as_match=treat_unknown_as_match)

    # TODO should have other functions, like printing at least...
    # MAYBE-TODO add overall info like total_reads, total_aligned, total_unaligned, total_groups, etc... 




### OLD STUFF FROM OLD VERSION, TODO REWRITE OR DELETE


def main(infiles, reffile, outfile, mutation_cutoffs, count_repeat_matches, allowed_repeats_file=None, _verbose=False):
    """ Read arguments, extract data, print output  - for details see module docstring."""
    global verbose
    verbose = _verbose
    if count_repeat_matches=="from_file" and not allowed_repeats_file:
        raise ValueError("If count_repeat_matches is from_file, an allowed_repeats_file argument is required.")
    if count_repeat_matches!="from_file" and allowed_repeats_file:
        print "Warning: The allowed_repeats_file will be ignored unless you set count_repeat_matches to from_file."
        allowed_repeats_file = None
    if count_repeat_matches not in ["all","none"] and not MATCH_TYPE_FIELD:
        print "Warning: you specified a count_repeat_matches method, but the MATCH_TYPE_FIELD variable in the program is set to 0, implying we're dealing with an input file format that doesn't provide match type information. All matches in the input file will be counted."
    if verbose: print "parsing reference file - time %s."%time.ctime()
    Name_to_seq, Name_to_mutations = initialize_dictionaries(reffile)
    if verbose: print "finished parsing reference file - time %s."%time.ctime()
    line_counts_summary = extract_data_from_novoalign_file(infiles, Name_to_seq, Name_to_mutations,
                                                           count_repeat_matches, allowed_repeats_file)
    total_lines = sum(line_counts_summary)
    output_summary = "### Finished! %i lines processed:\n"%total_lines
    output_summary += "# %i non-aligned lines, %i lines with non-allowed alignment types, %i alignments to reference file, %i alignments not matching reference file"%line_counts_summary
    if verbose: print "printing output - time %s."%time.ctime()
    Name_to_stats = calculate_statistics(Name_to_mutations)
    write_output(Name_to_seq,Name_to_stats,outfile,output_summary)
    if verbose: print output_summary

if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    # optparse only used here!  The options are converted to normal arguments before passing to main.
    from optparse import OptionParser
    USAGE = "USAGE: %prog [options] infile [more infiles] reference_file outfile"
    parser = OptionParser(usage = USAGE)
    parser.add_option('-m', '--mutation_cutoffs', default="1,3,10", metavar="<comma-separated-int-list>")
    parser.add_option('-r', '--count_repeat_matches', choices=["none","all","by_seq","from_file"], default="none", metavar="[none|all|by_seq|from_file]", help="How to deal with cases of one read aligning to multiple reference sequences: 'all' - count it as matching both sequences, 'none' - ignore it, 'by_seq' - only count it if both reference sequences are identical, 'from_file' - only count it if the pair occurs in the allowed_repeats_file (see -a option).")
    parser.add_option('-a', '--allowed_repeats_file', type="string", metavar="FILE", help="See -r option 'from_file' for the explanation. The file should contain any number of lines of the format \"sequence<tab> name1,name2,name3,...\", and may be generated by the find_library_duplicates.py script.")
    # TODO implement this with treating each mutated version of each shRNA as separate!
    parser.add_option('-M', '--separate_mutations', action="store_true", default=False, help="NOT IMPLEMENTED YET")
    parser.add_option('-v', '--verbose', action="store_true", default=False)
    (options, args) = parser.parse_args()
    # make sure required arguments are there
    if not len(args)>=3:    parser_error(parser,"infile, reference file and outfile are required")
    # parse mutation cutoff list
    try:                mutation_cutoffs = [int(x) for x in options.mutation_cutoffs.split(',')]
    except ValueError:  parser_error(parser,"mutation_cutoffs must be a comma-separated list of integers")
    # run the main program with the arguments as given (any number of infiles, one reffile, one outfile
    reffile, outfile = args[-2:]
    if len(args)==3:    infiles = [args[0]]
    else:               infiles = args[:-2]
    main(infiles, reffile, outfile, mutation_cutoffs, options.count_repeat_matches, options.allowed_repeats_file, options.verbose)

