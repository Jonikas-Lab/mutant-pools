#!/usr/bin/env python
"""
Module containing functions and classes useful for parsing deepseq SAM files and other deepseq data analysis we're doing.

This is a module to be imported and used by other programs.  Running it directly runs the built-in test suite.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011
"""

import sys
import unittest
from collections import defaultdict
from general_utilities import keybased_defaultdict
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
# MAYBE-TODO HTSeq doesn't appear aware of the = and X operations...  http://www-huber.embl.de/users/anders/HTSeq/doc/alignments.html#HTSeq.CigarOperation  - I emailed the author about it

VALID_POSITION_TYPES = ['leftmost','rightmost','5prime','3prime']


### Various functions

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
        # (need this []+ here so it's a copy, not a reference, and modifying it later doesn't modify the original)
        cigar_types_mutation = [] + CIGAR_TYPES_MUTATION
    else:
        cigar_types_mutation = CIGAR_TYPES_MUTATION + CIGAR_TYPES_INTRON   
    # figure out how to treat unknown matches ('M'), based on argument
    if treat_unknown_as=='unknown':
        # (need this []+ here so it's a copy, not a reference, and modifying it later doesn't modify the original)
        cigar_types_unknown = [] + CIGAR_TYPES_UNKNOWN
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
    #       basically a number means matches, a letter means a mismatch to reference (or insertion? is that different?), 
    #       letters preceded by ^ mean deletion from the reference
    try:                mutation_string = HTSeq_alignment.optional_field('MD')
    except KeyError:    return -1
    # for unalign reads MD field is missing - returns -1
    mutation_letters = [c for c in mutation_string if not (c.isdigit() or c=='^')]
    #   (^ is used in describing a mutation but it shouldn't be counted as a separate mutation - only letters count.)
    return len(mutation_letters)


def check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as='unknown', ignore_introns=False):
    """ Return number of mutations in HTSeq_alignment (look at CIGAR string and NM and MD optional fields); -1 if unknown.
    First check the CIGAR string but only accept the answer if there are no unknown ('M') characters; 
     then check the NM and MD fields and return the result if those fields exist.
    If the CIGAR string is ambiguous and neither of the optional fields exist:
     - if treat_unknown_as is 'unknown', return -1
     - if treat_unknown_as is 'mutation' or 'match', return the CIGAR string result with unknowns counted accordingly.
    If ignore_introns is False, count introns (N) in CIGAR string as mutations; otherwise don't .
    Does NOT guarantee returning a sensible value if the CIGAR, NM and MD fields contain inconsistent information.
    """
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
        raise ValueError("Invalid position %s! Need an HTSeq iv object. (If empty, maybe read wasn't aligned?)"%HTSeq_pos)
    chrom = HTSeq_pos.chrom
    if pos_type=='leftmost':        pos = HTSeq_pos.start
    elif pos_type=='rightmost':     pos = HTSeq_pos.end
    elif pos_type=='5prime':        pos = HTSeq_pos.start if HTSeq_pos.strand=='+' else HTSeq_pos.end
    elif pos_type=='3prime':        pos = HTSeq_pos.end if HTSeq_pos.strand=='+' else HTSeq_pos.start
    else:                           raise ValueError("pos_type argument must be one of %s."%VALID_POSITION_TYPES)
    return chrom, pos


### Main two classes

class Alignment_position_sequence_group():
    """ Data regarding sequences aligned to a particular genomic position (genomic position is set at initialization). 
    Variables: chromosome, position, total_read_count, perfect_read_count, unique_sequence_count, 
      sequences_and_counts (a sequence:count dictionary).
    Methods: add_read to add a given HTSeq read to the counts (doesn't check chromosome/position), 
     get_main_sequence to get the most common sequence from sequences_and_counts. """

    def __init__(self, chromosome, position):
        """ Set chromosome and position; initialize total_read_count, perfect_read_count and sequences_and_counts."""
        # should this check that chromosome is a string and position is an int?  
        #   Probably not, in case I want to use HTSeq or Biopython position objects later...
        self.chromosome, self.position = chromosome, position
        self.gene = 'gene_unknown'
        self.total_read_count = 0
        self.perfect_read_count = 0
        self.unique_sequence_count = 0
        self.sequences_and_counts = defaultdict(lambda: 0)

    # MAYBE-TODO give each mutant some kind of unique ID at some point in the process?  If we end up using per-mutant barcodes (in addition to the flanking sequences), we could use that, probably, or that plus genomic location

    def add_read(self, HTSeq_alignment, treat_unknown_as_match=False):
        """ Add a read to the data: increment total_read_count, increment perfect_read_count if read is a perfect 
        alignment, increment the appropriate field of sequences_and_counts based on read sequence."""
        # if it's a new sequence, increment unique_sequence_count; add a count to the self.sequences_and_counts dictionary.
        seq = HTSeq_alignment.read.seq
        if seq not in self.sequences_and_counts:
            self.unique_sequence_count += 1
        self.sequences_and_counts[seq] += 1
        # increment self.total_read_count; figure out if the read is perfect and increment self.perfect_read_count if yes.
        self.total_read_count += 1
        treat_unknown_as = 'match' if treat_unknown_as_match else 'mutation'
        mutation_count = check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as=treat_unknown_as)
        if mutation_count==0:  
            self.perfect_read_count += 1

    def add_counts(self, total_count, perfect_count, sequence_variant_count, assume_new_sequences=False):
        """ Increment self.total_read_count, self.perfect_read_count and self.unique_sequence_count based on inputs.
        Note that ifself.unique_sequence_count>0, it's impossible to determine the correct new value: 
         if we had old data with one unique sequence and now we have new data with another one, how do we know
          if that's the same or different sequence?  The correct total could be 1 or 2, so it's an option:
         If assume_new_sequences is True, the total is old+new; if it's False, the total is max(old,new).
        """
        self.total_read_count += total_count
        self.perfect_read_count += perfect_count
        if assume_new_sequences:
            self.unique_sequence_count += sequence_variant_count
        else:
            self.unique_sequence_count = max(self.unique_sequence_count,sequence_variant_count)

    def add_sequence_and_counts(self, seq, seq_count, add_to_uniqseqcount=True):
        """ Add seq_count to self.sequences_and_counts[seq] (it's created with count 0 if seq wasn't a key before).
        Note: if add_to_uniqseqcount is False, this will never increment self.unique_sequence_count;
         otherwise it only does so if seq was not already present in the self.sequences_and_counts data
          and if the total number of sequences in self.sequences_and_counts is higher than self.unique_sequence_count. """
        if add_to_uniqseqcount:
            if seq not in self.sequences_and_counts and len(self.sequences_and_counts)>self.unique_sequence_count:
                self.unique_sequence_count += 1
        self.sequences_and_counts[seq] += seq_count

    def get_main_sequence(self, N=1):
        """ Return the most common sequence in this group and its count (or Nth most common sequence if N is provided)."""
        sequences_by_count = sorted([(count,seq) for (seq,count) in self.sequences_and_counts.iteritems()], reverse=True)
        # try returning the Nth sequence and count; return nothing if there are under N sequences.
        try:                return (sequences_by_count[N-1][1], sequences_by_count[N-1][0])
        except IndexError:  return ('',0)


class All_alignments_grouped_by_pos():
    """ Essentially a dictionary of alignment_position_sequence_group with position data (chrom,pos) as keys. """

    def __init__(self, position_type=None):
        """ Checks position_type and assigns to self.position_type; initializes self.alignment_position_data_dict.
        position_type must be one of %s.
        self.alignment_position_data_dict is a dictionary that generates a new alignment_position_sequence_group object
         based on the key (i.e. position) if the key isn't already in the dictionary.  """%VALID_POSITION_TYPES
        self.alignment_position_data_dict = keybased_defaultdict(lambda key: Alignment_position_sequence_group(*key))
        if not position_type in VALID_POSITION_TYPES+[None]: 
            raise ValueError("The position_type variable must be one of %s!"%VALID_POSITION_TYPES)
        self.position_type = position_type
        self.total_read_count, self.aligned_read_count, self.unaligned_read_count = 0,0,0

    def add_alignment_reader_to_data(self, HTSeq_alignment_reader, treat_unknown_as_match=False):
        """ Adds all alignments to self.alignment_position_data_dict based on position.
        Input must be a list/generator/etc of HTSeq.Alignment objects (usually an HTSeq.SAM_Reader)."""
        if self.position_type is None:
            raise Exception("Cannot add data from an alignment reader if position_type isn't specified! Please set the position_type attribute of this All_alignments_grouped_by_pos instance to one of %s first."%VALID_POSITION_TYPES)
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
        position_info = "looking at %s end of read"%self.position_type
        OUTPUT.write("%s Read groups by alignment position (%s): %s\n"%(line_prefix, position_info, 
                                                                        len(self.alignment_position_data_dict)))

    def print_data(self, OUTPUT=None, sort_data=False, N_sequences=2, header_line=True, header_prefix="# "):
        """ Print the full data:  the read count for each group of sequences with the same position.
        If N_sequences<1, only prints position (chromosome and position), total and perfect read count, 
          and the number of sequence variants.  If N_sequences==1, also prints the most common sequence and count; 
          if N_sequences>1, adds the second most common sequence and count, and so on.
        Output is tab-separated, with headers starting with "# ".  Prints to an open file object (stdout by default).
        """
        if OUTPUT is None:
            OUTPUT = sys.stdout

        if header_line:
            headers = ['chromosome','position','gene','total_reads','perfect_reads', 'N_sequence_variants']
            for N in range(1,N_sequences+1):
                headers.extend(['sequence_%s'%N,'count_seq_%s'%N])
            OUTPUT.write(header_prefix + '\t'.join(headers) + "\n")

        if sort_data:
            data = sorted(self.alignment_position_data_dict.values(), key = lambda x: (x.chromosome, x.position))
        else:
            data = self.alignment_position_data_dict.itervalues()

        for group in data:
            group_data = [group.chromosome, group.position, group.gene, group.total_read_count, 
                          group.perfect_read_count, len(group.sequences_and_counts)]
            OUTPUT.write('\t'.join([str(x) for x in group_data]))
            for N in range(1,N_sequences+1):
                OUTPUT.write('\t%s\t%s'%group.get_main_sequence(N))
                # MAYBE-TODO also give the length and number of mutations for each sequence? Optionally?  Length is easy, but do I even keep track of mutation number?  I probably should...
            OUTPUT.write('\n')

    def read_from_file(self, infile, assume_new_sequences=False):
        """ Read data from a file generated by self.print_data, and add to self.alignment_position_data_dict.  
        Note that if you read data from a file and add it to data that already exists, it's impossible to determine the 
         correct number of unique sequence variants per group: if the data originally listed 1 unique sequence, and new
          data read from a file adds another 2 sequences, is that a total of 2 or 3 unique sequences? 
         If assume_new_sequences is True, the total is old+new; if it's False, the total is max(old,new).
        Currently doesn't read the specific sequences and counts - even if it did, some information could 
        always be missing, since the file only has N first sequences. """
        for line in open(infile):
            # MAYBE-TODO get unaligned read count from summary if present - do we actually want that?
            # ignore comment and header lines
            if line.startswith('#'):    continue
            if line.startswith('chromosome\tposition\t'):    continue
            fields = line.split('\t')
            chromosome = fields[0]
            pos = int(fields[1])
            gene = fields[2]
            total_reads,perfect_reads,sequence_variants = [int(x) for x in fields[3:6]]
            self.alignment_position_data_dict[(chromosome,pos)].add_counts(total_reads, perfect_reads, 
                                                                           sequence_variants, assume_new_sequences)
            self.alignment_position_data_dict[(chromosome,pos)].gene = gene
            # MAYBE-TODO split the remaining fields into twos (seq,count) (we don't know how many there will be) 
            #   and use self.alignment_position_data_dict[(chromosome, position)].add_sequence_and_counts(seq,count) 
            #   to add specific sequences and counts to the data.  Add to unit-test!
        self.aligned_read_count = sum([x.total_read_count for x in self.alignment_position_data_dict.values()])
        self.total_read_count = self.aligned_read_count + self.unaligned_read_count
    

# MAYBE-TODO at some point I'll probably want to generalize both of those classes to not necessarily be grouped by position...


######### Test code #########

class Testing_single_functions(unittest.TestCase):
    """ Unit-tests for all the top-level functions in this module. """

    class Fake_cigar_op:
        """ Fake CIGAR operation, mimicking HTSeq cigar object."""
        def __init__(self,string):  self.type = string
        size = 1
    class Fake_alignment:
        """ Fake alignment with optional_field lookup function, mimicking HTSeq alignment object."""
        def __init__(self,data):        self.data = data
        def optional_field(self,x):     return self.data[x]

    def test__check_mutation_count_by_CIGAR_string(self):
        # no alignment (CIGAR is None)
        class fake_alignment:
            cigar = None
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == -1
        # CIGAR is unambiguous, no MD or NM given (or needed)
        class fake_alignment:
            cigar = [self.Fake_cigar_op('='),self.Fake_cigar_op('=')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 0
        class fake_alignment:
            cigar = [self.Fake_cigar_op('X'),self.Fake_cigar_op('X')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 2
        class fake_alignment:
            cigar = [self.Fake_cigar_op('D'),self.Fake_cigar_op('D')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 2
        class fake_alignment:
            cigar = [self.Fake_cigar_op('S'),self.Fake_cigar_op('=')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 0
        class fake_alignment:
            cigar = [self.Fake_cigar_op('S'),self.Fake_cigar_op('X')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 1
        class fake_alignment:
            cigar = [self.Fake_cigar_op('N'),self.Fake_cigar_op('=')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment, ignore_introns=True) == 0
        assert check_mutation_count_by_CIGAR_string(fake_alignment, ignore_introns=False) == 1
        # CIGAR is ambiguous (contains M's) - return -1, 2 or 0 depending on what treat_unknown_as is set to
        class fake_alignment:
            cigar = [self.Fake_cigar_op('M'),self.Fake_cigar_op('M')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 2
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 0
        class fake_alignment:
            cigar = [self.Fake_cigar_op('M'),self.Fake_cigar_op('=')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 0
        class fake_alignment:
            cigar = [self.Fake_cigar_op('M'),self.Fake_cigar_op('X')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 2
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 1

    def test__check_mutation_count_by_optional_NM_field(self):
        """ the tested function should return -1 if no NM field, otherwise return value of NM field. """
        fake_alignment = self.Fake_alignment({})
        assert check_mutation_count_by_optional_NM_field(fake_alignment) == -1
        for x in range(10):
            fake_alignment = self.Fake_alignment({'NM':x})
            assert check_mutation_count_by_optional_NM_field(fake_alignment) == x

    def test__check_mutation_count_by_optional_MD_field(self):
        """ see ~/experiments/reference_data/aligner_format_info/* files for MD field examples."""
        fake_alignment = self.Fake_alignment({})
        assert check_mutation_count_by_optional_MD_field(fake_alignment) == -1
        for s in [str(x) for x in range(30)]:
            fake_alignment = self.Fake_alignment({'MD': s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 0
            fake_alignment = self.Fake_alignment({'MD': s+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 0
            fake_alignment = self.Fake_alignment({'MD': s+'A'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 1
            fake_alignment = self.Fake_alignment({'MD': s+'A0G'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 2
            fake_alignment = self.Fake_alignment({'MD': s+'A2G'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 2
            fake_alignment = self.Fake_alignment({'MD': s+'A2G2T2C2N'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 5
            fake_alignment = self.Fake_alignment({'MD': s+'^A'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 1
            fake_alignment = self.Fake_alignment({'MD': s+'^AGC'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 3

    def test__check_mutation_count_try_all_methods(self):
        """ The order of check is CIGAR, NM, MD; CIGAR is skipped if ambiguous; NM and MD skipped if inexistent. 
        Not attempting to deal with inconsistent states sensibly."""
        # all measures agree there are no mutations (with 0-2 of NM/MD fields present)
        for optional_data in [{'NM':0, 'MD':'10'}, {'NM':0}, {'MD':'10'}, {}]:
            fake_alignment = self.Fake_alignment(optional_data)
            fake_alignment.cigar = [self.Fake_cigar_op('=') for x in range(10)]
            assert check_mutation_count_try_all_methods(fake_alignment) == 0
        # all measures agree there is a mutation (with 0-2 of NM/MD fields present)
        for optional_data in [{'NM':1, 'MD':'A9'}, {'NM':1}, {'MD':'A9'}, {}]:
            fake_alignment = self.Fake_alignment(optional_data)
            fake_alignment.cigar = [self.Fake_cigar_op('X')] + [self.Fake_cigar_op('=') for x in range(9)]
            assert check_mutation_count_try_all_methods(fake_alignment) == 1
        # CIGAR is ambiguous, there are no mutations according to NM/MD (NM, MD or both are present)
        for optional_data in [{'NM':0, 'MD':'10'}, {'NM':0}, {'MD':'10'}]:
            fake_alignment = self.Fake_alignment(optional_data)
            fake_alignment.cigar = [self.Fake_cigar_op('M') for x in range(10)]
            assert check_mutation_count_try_all_methods(fake_alignment) == 0
        # CIGAR is ambiguous, there is a  mutation according to NM/MD (NM, MD or both are present)
        for optional_data in [{'NM':1, 'MD':'A9'}, {'NM':1}, {'MD':'A9'}]:
            fake_alignment = self.Fake_alignment(optional_data)
            fake_alignment.cigar = [self.Fake_cigar_op('M') for x in range(10)]
            assert check_mutation_count_try_all_methods(fake_alignment) == 1

    class Fake_pos:
        def __init__(self, chrom, strand, start, end):
            self.chrom = chrom
            self.strand = strand
            self.start = start
            self.end = end

    def test__get_chrom_and_pos_from_HTSeq_position(self):
        """ Expected results: leftmost=start, rightmost=end; 5prime=start and 3prime=end if strand is +, reverse otherwise.
        Tested function should raise ValueError when passed None or an invalid pos_type. """
        for pos_type in VALID_POSITION_TYPES:
            self.assertRaises(ValueError, get_chrom_and_pos_from_HTSeq_position, None, pos_type)
        fake_pos = self.Fake_pos('C', '+', 0, 5)
        for bad_pos_type in ['', 'aaa', 0, 1, [], None, True, False, 'start', 'end', 'middle', 'read']:
            self.assertRaises(ValueError, get_chrom_and_pos_from_HTSeq_position, fake_pos, bad_pos_type)
        for (start,end) in [(0,5), (0,100), (10,11), (5,44)]:
            fake_pos = self.Fake_pos('C', '+', start, end)
            assert get_chrom_and_pos_from_HTSeq_position(fake_pos, 'leftmost') == ('C', start)
            assert get_chrom_and_pos_from_HTSeq_position(fake_pos, 'rightmost') == ('C', end)
            assert get_chrom_and_pos_from_HTSeq_position(fake_pos, '5prime') == ('C', start)
            assert get_chrom_and_pos_from_HTSeq_position(fake_pos, '3prime') == ('C', end)
            fake_pos = self.Fake_pos('C', '-', start, end)
            assert get_chrom_and_pos_from_HTSeq_position(fake_pos, 'leftmost') == ('C', start)
            assert get_chrom_and_pos_from_HTSeq_position(fake_pos, 'rightmost') == ('C', end)
            assert get_chrom_and_pos_from_HTSeq_position(fake_pos, '5prime') == ('C', end)
            assert get_chrom_and_pos_from_HTSeq_position(fake_pos, '3prime') == ('C', start)


class Testing_Alignment_position_sequence_group(unittest.TestCase):
    """ Unit-tests for the Alignment_position_sequence_group class and its methods. """

    def test__init(self):
        for chromosome in ['chr1', 'chromosome_2', 'chrom3', 'a', 'adfads', '100', 'scaffold_88']:
            for position in [0,1,2,5,100,10000,4323423]:
                group = Alignment_position_sequence_group(chromosome,position)
                assert group.chromosome == chromosome
                assert group.position == position
                assert group.total_read_count == 0
                assert group.perfect_read_count == 0
                assert group.unique_sequence_count == 0
                assert group.sequences_and_counts == {}

    def test__add_read(self):
        pass
        # MAYBE-TODO implement using a mock-up of HTSeq_alignment?  (see Testing_single_functions for how I did that)

    def test__add_counts(self):
        group = Alignment_position_sequence_group('chr',3)
        group.add_counts(0,0,0)
        assert group.total_read_count == 0
        assert group.perfect_read_count == 0
        assert group.sequences_and_counts == {}
        group.add_counts(2,2,1)
        assert group.total_read_count == 2
        assert group.perfect_read_count == 2
        assert group.sequences_and_counts == {}
        group.add_counts(2,1,1)
        assert group.total_read_count == 4
        assert group.perfect_read_count == 3
        assert group.sequences_and_counts == {}
        # how group.unique_sequence_count changes depends on assume_new_sequences:
        #  - if False, group.unique_sequence_count is max(new_seq_count, group.unique_sequence_count)
        group = Alignment_position_sequence_group('chr',3)
        group.add_counts(0,0,0,assume_new_sequences=False)
        assert group.unique_sequence_count == 0
        group.add_counts(1,1,1,assume_new_sequences=False)
        assert group.unique_sequence_count == 1
        group.add_counts(2,2,2,assume_new_sequences=False)
        assert group.unique_sequence_count == 2
        group.add_counts(1,1,1,assume_new_sequences=False)
        assert group.unique_sequence_count == 2
        group.add_counts(2,2,2,assume_new_sequences=False)
        assert group.unique_sequence_count == 2
        #  - if True and group.unique_sequence_count is new_seq_count + group.unique_sequence_count
        group = Alignment_position_sequence_group('chr',3)
        group.add_counts(0,0,0,assume_new_sequences=True)
        assert group.unique_sequence_count == 0
        group.add_counts(1,1,1,assume_new_sequences=True)
        assert group.unique_sequence_count == 1
        group.add_counts(2,2,2,assume_new_sequences=True)
        assert group.unique_sequence_count == 3
        group.add_counts(1,1,1,assume_new_sequences=True)
        assert group.unique_sequence_count == 4
        group.add_counts(2,2,2,assume_new_sequences=True)
        assert group.unique_sequence_count == 6
        group.add_counts(2,2,2,assume_new_sequences=False)
        assert group.unique_sequence_count == 6

    def test__add_sequence_and_counts(self):
        group = Alignment_position_sequence_group('chr',3)
        # adding sequence/count to group.sequences_and_counts, WITHOUT touching group.unique_sequence_count
        group.add_sequence_and_counts('AAA',2,add_to_uniqseqcount=False)
        assert group.sequences_and_counts == {'AAA':2}
        assert group.unique_sequence_count == 0
        # adding sequence/count to group.sequences_and_counts, and incrementing group.unique_sequence_count if warranted:
        #  - if adding a sequence that was already there, don't increment
        group.add_sequence_and_counts('AAA',2,add_to_uniqseqcount=True)
        assert group.sequences_and_counts == {'AAA':4}
        assert group.unique_sequence_count == 0
        #  - if adding a new sequence, increment
        group.add_sequence_and_counts('GGG',2,add_to_uniqseqcount=True)
        assert group.sequences_and_counts == {'AAA':4, 'GGG':2}
        assert group.unique_sequence_count == 1
        #  - if adding a new sequence but group.unique_sequence_count is already higher than expected, don't increment
        group.unique_sequence_count = 5
        group.add_sequence_and_counts('CCC',2,add_to_uniqseqcount=True)
        assert group.sequences_and_counts == {'AAA':4, 'GGG':2, 'CCC':2}
        assert group.unique_sequence_count == 5
        # make sure it raises an error if given a non-numeric argument
        for not_a_number in [None,'','a','GGG']:
            self.assertRaises(TypeError,group.add_sequence_and_counts,'CCC',not_a_number)

    def test__get_main_sequence(self):
        group = Alignment_position_sequence_group('chr',3)
        assert group.get_main_sequence() == ('',0)
        assert group.get_main_sequence(1) == ('',0)
        assert group.get_main_sequence(4) == ('',0)
        group.add_sequence_and_counts('AAA',1)
        group.add_sequence_and_counts('GGG',2)
        assert group.sequences_and_counts == {'AAA':1, 'GGG':2}
        assert group.get_main_sequence() == ('GGG',2)
        assert group.get_main_sequence(1) == ('GGG',2)
        assert group.get_main_sequence(2) == ('AAA',1)
        assert group.get_main_sequence(3) == ('',0)
        assert group.get_main_sequence(4) == ('',0)
        group.add_sequence_and_counts('CCC',1)
        group.add_sequence_and_counts('AAA',2)
        assert group.sequences_and_counts == {'AAA':3, 'GGG':2, 'CCC':1}
        assert group.get_main_sequence() == ('AAA',3)
        assert group.get_main_sequence(1) == ('AAA',3)
        assert group.get_main_sequence(2) == ('GGG',2)
        assert group.get_main_sequence(3) == ('CCC',1)
        assert group.get_main_sequence(4) == ('',0)
        assert group.get_main_sequence(5) == ('',0)


class Testing_All_alignments_grouped_by_pos(unittest.TestCase):
    """ Unit-tests for the All_alignments_grouped_by_pos class and its methods. """

    def test__init(self):
        for pos_type in VALID_POSITION_TYPES+[None]:
            data = All_alignments_grouped_by_pos(pos_type)
            assert data.position_type == pos_type
            assert data.alignment_position_data_dict == {}
            assert data.total_read_count == 0
            assert data.aligned_read_count == 0
            assert data.unaligned_read_count == 0
        for pos_type in [True, False, 0, 0.11, 23, 'asdfas', '', 'something', [2,1]]:
            self.assertRaises(ValueError, All_alignments_grouped_by_pos, pos_type)

    def test__add_alignment_reader_to_data(self):
        pass
        # TODO make sure it fails if self.position_type isn't defined...
        # MAYBE-TODO implement using a mock-up of HTSeq_alignment?  (see Testing_single_functions for how I did that)

    def test__print_summary(self):
        pass
        # MAYBE-TODO implement based on stuff in test_data, like do_test_run in deepseq_count_alignments.py?

    def test__print_data(self):
        pass
        # MAYBE-TODO implement based on stuff in test_data, like do_test_run in deepseq_count_alignments.py?

    def test__read_from_file(self):
        input_file = 'test_data/test_output__U_leftmost.txt'
        data = All_alignments_grouped_by_pos(None)
        data.read_from_file(input_file)
        assert data.aligned_read_count == 29
        assert data.unaligned_read_count == 0
        assert data.total_read_count == 29
        # just spot-checking some of the outputs
        group = data.alignment_position_data_dict[('reads_2_seqs_1',199)]
        assert group.chromosome == 'reads_2_seqs_1'
        assert group.position == 199
        assert group.total_read_count == 2
        assert group.perfect_read_count == 2 
        assert group.unique_sequence_count == 1
        group = data.alignment_position_data_dict[('mutation_yes',199)]
        assert group.chromosome == 'mutation_yes'
        assert group.position == 199
        assert group.total_read_count == 6
        assert group.perfect_read_count == 0
        assert group.unique_sequence_count == 1
        group = data.alignment_position_data_dict[('strandedness_+_reference',99)]
        assert group.chromosome == 'strandedness_+_reference'
        assert group.position == 99
        assert group.total_read_count == 1
        assert group.perfect_read_count == 1
        assert group.unique_sequence_count == 1
        # try adding more data to a file that already has some...
        data.read_from_file(input_file, assume_new_sequences=False)
        group = data.alignment_position_data_dict[('reads_2_seqs_1',199)]
        assert group.chromosome == 'reads_2_seqs_1'
        assert group.position == 199
        assert group.total_read_count == 4
        assert group.perfect_read_count == 4
        # how group.unique_sequence_count should act in this case depends on the value of assume_new_sequences
        assert group.unique_sequence_count == 1
        data.read_from_file(input_file, assume_new_sequences=True)
        assert group.unique_sequence_count == 2


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main(argv=[sys.argv[0]])
    # if I wanted more control I could do this instead:
    #import os
    #unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromName(os.path.splitext(sys.argv[0])[0]))
    #   (autodetection of all tests - see http://docs.python.org/library/unittest.html#unittest.TestLoader)
    # there's probably also some way to easily get all tests from the current file without passing the name, but I haven't found it yet...

