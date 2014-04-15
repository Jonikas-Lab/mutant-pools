#!/usr/bin/env python2.7
"""
Module containing classes and functions for analysis of deepseq data related to insertional mutant libraries.

This is a module to be imported and used by other programs.  Running it directly runs the built-in unit-test suite.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2013
"""

from __future__ import division
# basic libraries
import sys, re
import unittest
from collections import defaultdict, Counter
from math import isnan
# other libraries
import HTSeq
import matplotlib.pyplot as mplt
import plotting_utilities
# my modules
from mutant_analysis_classes import MutantError, Insertion_position, Insertional_mutant, Insertional_mutant_pool_dataset, is_cassette_chromosome, HTSeq_pos_to_tuple, SEQ_ENDS, SEQ_STRANDS
from general_utilities import keybased_defaultdict, value_and_percentages
from seq_basic_utilities import parse_fasta
from deepseq_utilities import Fake_deepseq_objects
# there's a "from parse_annotation_file import parse_gene_annotation_file" in one function, not always needed


# TODO need to implement new Insertion_position variant with uncertainty!  And modify the gene/feature finding function to deal with that...

############################### Constants and help functions #####################################

MAX_POSITION_DISTANCE = 1500

def get_Carette_pos_from_flanking_region_pos(flanking_region_pos, cassette_end, reads_are_reverse=False, immutable_position=True):
    """ Return a Insertion_position instance giving the position of the far end of Carette genome-side read.

    The strand should be the same as that for the cassette-side read if the two are consistent.

    Flanking_region_Pos should be a HTSeq.GenomicPosition instance, or a (chrom,start_pos,end_pos,strand) tuple
      (the tuple should have 1-based end-inclusive positions, so AA is 1-2 in AATT; HTSeq positions are 0-based end-exclusive); 
     cassette_end gives the side of the insertion that read is on; 
     reads_are_reverse is True if the read is in reverse orientation to the cassette, False otherwise. 

    See mutant_analysis_classes.get_insertion_pos_from_flanking_region_pos for how this is done for insertion positions.
    This is somewhat different: we know the genome-side read is in opposite orientation from the cassette-side read
     (since they're paired-end reads from the same data), so in order for the strand to match between these two ends, 
     it has to be calculated in the opposite way. 
    The position we're interested in is always the position of the START of the read 
     (which according to most encodings is the end if the strand is -).

    If immutable_position is True, the position will be made immutable after creation (this is reversible).
    """
    # note: For Carette cassette-side reads, 5' data has 5prime and reads_are_reverse=True, 3' data has 3prime and reads_are_reverse=False.  NOT DEALING WITH other cases!
    try:                                chrom, start_pos, end_pos, strand = flanking_region_pos
    except (TypeError, ValueError):     chrom, start_pos, end_pos, strand = HTSeq_pos_to_tuple(flanking_region_pos) 
    if strand not in SEQ_STRANDS:       raise MutantError("Invalid strand %s! Must be one of %s"%(strand, SEQ_STRANDS))
    if cassette_end not in SEQ_ENDS:    raise MutantError("Invalid end %s! Must be one of %s"%(cassette_end, SEQ_ENDS))
    if start_pos < 1:                   raise MutantError("Flanking region positions must be positive!")
    if start_pos > end_pos:             raise MutantError("Flanking region start can't be after end!")
    if not ((cassette_end=='5prime' and reads_are_reverse==True) or (cassette_end=='3prime' and reads_are_reverse==False)):
        raise MutantError("For a Carette read, reads have to read OUT of the cassette - check your end/direction!")
    ### chromosome is always the same as read, so just leave it as is
    ### we always want min_position to be the position of the start of the read, 
    #     since that's the "outer" position in a paired-end set.
    if strand=='+':     pos_before, pos_after = start_pos, None
    else:               pos_before, pos_after = end_pos, None
    ### cassette strand is the opposite to how it is with the cassette-side read
    if cassette_end=='3prime':  strand = '-' if strand=='+' else '+'
    return Insertion_position(chrom, strand, position_before=pos_before, position_after=pos_after, immutable=immutable_position)


############################### Main classes describing the mutants and mutant sets #####################################

# help functions returning "blank" mutants for defaultdicts (can't use lambdas because pickle doesn't like them)
def blank_Carette_mutant():
    return Insertional_mutant_Carette()
def blank_Carette_mutant_with_pos(pos):
    return Insertional_mutant_Carette(insertion_position=pos)


class Insertional_mutant_Carette(Insertional_mutant):
    # TODO docstring!
    # TODO unit-test all this!!

    def __init__(self, *args, **kwargs):
        """ Same as for Insertional_mutant, plus add empty set of read IDs. """
        Insertional_mutant.__init__(self, *args, **kwargs)
        self.read_IDs = set()
        self.Carette_genome_side_reads = []

    def add_read(self, HTSeq_alignment, read_count=1, treat_unknown_as_match=False, dataset_name=None):
        """ Add a read to the data (or multiple identical reads, if read_count>1); return True if perfect match.

        Runs Insertional_mutant.add_read, plus adds read ID to self.read_IDs.
        """
        Insertional_mutant.add_read(self, HTSeq_alignment, read_count, treat_unknown_as_match, dataset_name)
        self.read_IDs.add(HTSeq_alignment.read.name)

    def merge_mutant(self, *args, **kwargs):
        raise MutantError("merge_mutant NOT IMPLEMENTED on Carette mutant object!")
        # TODO implement?  Just need to merge the read IDs and Carette read data...

    def add_Carette_read(self, new_position, HTSeq_alignment=None, seq=None, N_errors=None, read_count=1):
        """ Must provide either HTSeq_alignment (if aligned) or seq (if not).
        If providing HTSeq_alignment, seq and N_errors will be ignored. 
        Must provide new_position (which should be an Insertion_position instance if read is uniquely aligned, 
         or None/"Multiple" or such if not); if HTSeq_alignment is given, method does NOT check that new_position matches it.
        """
        # TODO make sure this makes sense!
        mutant = Insertional_mutant(new_position)
        if HTSeq_alignment is not None and seq is not None:
            raise MutantError("When running add_Carette_read, can't give BOTH HTSeq_alignment and seq!")
        elif HTSeq_alignment is not None:
            mutant.add_read(HTSeq_alignment, read_count)
        elif seq is not None:
            mutant.add_sequence_and_counts(seq, read_count)
            if N_errors==0:
                mutant.add_counts(0, read_count, 0)
        else:
            raise MutantError("When running add_Carette_read, must give HTSeq_alignment or seq!")
        self.Carette_genome_side_reads.append(mutant)
        # Note: adding gene/annotation info for those is implemented in the dataset methods.

    def _if_confirming_read(self, position, max_distance):
        # read needs to match chromosome and strand to be consistent
        if position.chromosome != self.position.chromosome:                         return False
        if position.strand != self.position.strand:                                 return False
        # it also needs to be below max_distance
        if abs(position.min_position - self.position.min_position) > max_distance:  return False
        return True
        # TODO what about WEIRD cases, when it's the right distance but in the wrong direction? Those are rare enough to ignore for now, but should do something about them!

    def _decide_if_replace_read(self, new_position, max_distance):
        """ Decide whether new position is "better" than old one - just refactored out of improve_best_Carette_read for convenience.
        """
        # if there are no current reads, add new one
        if not len(self.Carette_genome_side_reads):                     return True
        # if new read isn't "confirming", it can't be better
        if not self._if_confirming_read(new_position, max_distance):    return False
        # if the new one is "confirming" and the old one isn't, new one has to be better
        old_position = self.Carette_genome_side_reads[0].position
        if not self._if_confirming_read(old_position, max_distance):    return True
        # if both the old and new position meet the basic conditions, pick the highest-distance one
        new_dist = abs(new_position.min_position - self.position.min_position)
        old_dist = abs(old_position.min_position - self.position.min_position)
        if new_dist > old_dist:                                 return True
        else:                                                   return False

    def improve_best_Carette_read(self, new_position, HTSeq_alignment, max_distance=MAX_POSITION_DISTANCE):
        """ Compare the read to the current best one - replace the best one if this is better.

        "Better" meaning furthest away from the cassette-side read, while still remaining on the same chromosome, strand,
            and within max_distance of it.
         If both reads are on a different chromosome or strand than the cassette-side position (both are bad, really), 
            the first one is kept.
         If there is no current read, the new one is always used.
        """
        # if there are more than one current reads, you're not using improve_best_Carette_read consistently!
        if len(self.Carette_genome_side_reads) > 1:
            raise MutantError("Don't try using the improve_best_Carette_read when keeping more than one read!")
        # if decided to replace, just make self.Carette_genome_side_reads be the new mutant, discarding old version
        if self._decide_if_replace_read(new_position, max_distance):
            new_mutant = Insertional_mutant(new_position)
            new_mutant.add_read(HTSeq_alignment)
            self.Carette_genome_side_reads = [new_mutant]

    def Carette_N_confirming_reads(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return the number of aligned genome-side reads that DON'T confirm the cassette-side position.

        ("Confirm" means same chromosome, consistent strand, and at most max_distance away)
        """
        N = 0
        for read_data in self.Carette_genome_side_reads:
            # skip non-aligned reads
            try:                    chrom = read_data.position.chromosome
            except AttributeError:  continue
            if self._if_confirming_read(read_data.position, max_distance):  N += 1
        return N

    def Carette_N_non_confirming_reads(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return the number of aligned genome-side reads that DON'T confirm the cassette-side position.

        ("Confirm" means same chromosome, consistent strand, and at most max_distance away)
        """
        return self.Carette_N_genomic_reads + self.Carette_N_cassette_reads - self.Carette_N_confirming_reads(max_distance)

    @property
    def Carette_N_unaligned_reads(self):
        return sum(1 for read_data in self.Carette_genome_side_reads if read_data.position in 'UNALIGNED MULTIPLE'.split())
        # TODO if I change to treating multiple as aligned, will this need to change?
        # TODO all this UNALIGNED/MULTIPLE stuff is text, should probably change it to constants!

    @property
    def Carette_N_genomic_reads(self):
        N = 0
        for read_data in self.Carette_genome_side_reads:
            try:                    chrom = read_data.position.chromosome
            except AttributeError:  continue
            if not is_cassette_chromosome(chrom):
                N += 1
        return N

    @property
    def Carette_N_cassette_reads(self):
        return len(self.Carette_genome_side_reads) - self.Carette_N_genomic_reads - self.Carette_N_unaligned_reads

    @property
    def Carette_N_genomic_chromosomes(self):
        chroms = set()
        for read_data in self.Carette_genome_side_reads:
            # grab chromosome - unless read is unaligned, then ignore and go on to next one
            try:                    chrom = read_data.position.chromosome
            except AttributeError:  continue
            chroms.add(chrom)
        return sum(1 for chrom in chroms if not is_cassette_chromosome(chrom))
            
    def Carette_N_distinct_positions(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return number of distinct insertion positions implied by Carette data (genomic, cassette, and #chromosomes).

        The output is a 3-tuple giving the number of distinct genome and cassette positions, 
         and the number of distinct non-cassette chromosomes the positions are in.
        Positions on different chromosomes/strands are always counted as distinct; positions on same chromosome/strand are 
         counted as distinct if the distance between them is >=max_distance (THIS IS SLIGHTLY ROUGH).
        
        Data used to generate positions includes the cassette-side position (single) and all the genome-side Carette read positions. 
        Unaligned reads are ignored.
        """
        positions_by_chrom_strand = defaultdict(list)
        # add the cassette-side position (single)
        positions_by_chrom_strand[(self.position.chromosome, self.position.strand)].append(self.position.min_position)
        # add all the genome-side read positions; skip unaligned ones.
        for read_data in self.Carette_genome_side_reads:
            pos = read_data.position
            try:                    positions_by_chrom_strand[(pos.chromosome, pos.strand)].append(pos.min_position)
            except AttributeError:  continue
        # count total number of dictinct positions - different chromosomes or strands, or distance > max_distance
        total_distinct_positions_genome = 0
        total_distinct_positions_cassette = 0
        # for each chromosome, go over all positions and only count ones every MAX_POSITION_DISTANCE as distinct
        for chrom_strand, positions in positions_by_chrom_strand.items():
            positions.sort()
            distinct_positions = [positions[0]]
            for pos in positions[1:]:
                if (pos-distinct_positions[-1]) > max_distance:
                    distinct_positions.append(pos)
            if is_cassette_chromosome(chrom_strand[0]):     total_distinct_positions_cassette += len(distinct_positions)
            else:                                           total_distinct_positions_genome += len(distinct_positions)
        return total_distinct_positions_genome, total_distinct_positions_cassette
    
    def Carette_confirmation_status(self, side=None, max_distance=MAX_POSITION_DISTANCE, min_weird_distance=20):
        """ Return status based on comparison of cassette-side and genome-side reads.

        Also checks whether results make sense given which side of the cassette we're looking at, if provided
        """
        # TODO better docstring, unit-tests!
        if not self.total_read_count:
            return 'NOT_FOUND'
        elif not self.Carette_N_genomic_reads + self.Carette_N_cassette_reads:
            return 'NO_ALIGNED_READS'
        # We take "sum(mutant.Carette_N_distinct_positions(max_allowed_distance)) > 1" to check for EITHER genomic OR cassette positions - TODO rewrite that method to be more obvious?
        elif sum(self.Carette_N_distinct_positions(max_distance)) > 1:
            return 'WRONG_POSITION'
        else:
            distances = []
            for Carette_read_data in self.Carette_genome_side_reads:
                try:                    distances.append(Carette_read_data.position.min_position - self.position.min_position)
                except AttributeError:  pass
            # check to make sure the position of the genomic-side reads MAKES SENSE vs the cassette-side 
            #  (is at least M earlier for +strand 5' and -strand 3', and M later for the other two cases)
            # TODO really, min_weird_distance default should probably depend on the read length...?
            if (((self.position.strand=='+' and side=="5'") or (self.position.strand=='-' and side=="3'")) 
                and max(distances) > -min_weird_distance):
                print ("WEIRD case: a genome-side Carette read starts only %sbp from the insertion site! %s")%(-max(distances), 
                                                                                                             self.position)
                return 'WEIRD'
            elif (((self.position.strand=='-' and side=="5'") or (self.position.strand=='+' and side=="3'")) 
                  and min(distances) < min_weird_distance):
                print ("WEIRD case: a genome-side Carette read starts only %sbp from the insertion site! %s")%(min(distances), 
                                                                                                             self.position)
                return 'WEIRD'
            else:
                return 'CONFIRMED'

    def Carette_max_confirmed_distance(self, max_distance=MAX_POSITION_DISTANCE):
        """ Return the distance between the cassette-side read and the furthest-away same-area genome-side read.

        Return NaN if no same-area genome-side reads.
        """
        # TODO better docstring, unit-tests!
        distances = []
        for Carette_read_data in self.Carette_genome_side_reads:
            # Only look at the genome-side reads that match the cassette-side read position!
            # There's a try/except because unaligned reads don't have proper positions.
            try:                    
                if (Carette_read_data.position.chromosome == self.position.chromosome 
                    and Carette_read_data.position.strand == self.position.strand):
                    pos_difference = abs(Carette_read_data.position.min_position - self.position.min_position)
                else:
                    continue
            except AttributeError:  
                continue
            if pos_difference <= max_distance:
                distances.append(pos_difference)
        try:
            return max(distances)
        except ValueError:
            return float('NaN')

    # TODO add new method that infers the approximate real insertion location!

    def Carette_print_detail(self, OUTPUT=sys.stdout, max_distance=MAX_POSITION_DISTANCE):
        ### print summary line
        N_distinct_genome, N_distinct_cassette = self.Carette_N_distinct_positions(max_distance)
        OUTPUT.write(" * %s distinct genomic positions (on %s chromosomes; %s genomic reads), "%(N_distinct_genome, 
                                                                 self.Carette_N_genomic_chromosomes, self.Carette_N_genomic_reads)
                     +"plus %s distinct cassette positions (%s reads) and %s unaligned/multi-aligned reads:\n"%(N_distinct_cassette,
                                                                 self.Carette_N_cassette_reads, self.Carette_N_unaligned_reads))
        # TODO is that a useful summary?  Is it confusing that the casette-side read isn't included in the read numbers?
        ### print cassette-side position
        OUTPUT.write("Cassette-side position:::\n")
        main_pos_fields = [self.position.chromosome, self.position.strand, self.position.full_position, self.gene, self.orientation, 
                           self.gene_feature, self.get_main_sequence()[0]] + self.gene_annotation
        OUTPUT.write('\t'.join([str(x) for x in main_pos_fields]) + '\n')
        ### print lines for each of the Carette genome-side reads
        # sort the Carette reads by alignment position (chrom,strand,pos; unaligned go last)
        OUTPUT.write("Carette protocol genome-side reads:::\n")
        for read_data in sorted(self.Carette_genome_side_reads, key = lambda x: x.position):
            try: 
                fields = [read_data.position.chromosome, read_data.position.strand, read_data.position.min_position]
            except AttributeError:
                fields = [read_data.position, '-', '-']
            fields += [read_data.gene, read_data.orientation, read_data.gene_feature, read_data.get_main_sequence()[0]] 
            try:
                fields += read_data.gene_annotation
            except AttributeError:
                pass
            OUTPUT.write('\t'.join([str(x) for x in fields]) + '\n')

    # TODO TODO TODO finish implementing class!  What other methods to add/overwrite?

# TODO what about mutants that combine normal MmeI and Carette data - are those going to be any different?


class Insertional_mutant_pool_dataset_Carette(Insertional_mutant_pool_dataset):
    # TODO docstring!

    # TODO this needs proper extra summary stuff! 
    #   #correct/incorrect positions, #cassette fragments, for correct positions distribution of distance checked, ...

    def __init__(self, cassette_end='?', reads_are_reverse='?', multi_dataset=False, infile=None):
        """ Initializes empty dataset; saves properties as provided; optionally reads data from mutant infile. """
        # Some things aren't implemented for Carette datasets
        if multi_dataset:
            raise MutantError("Multi-dataset Carette mutant pool datasets not implemented!")
        if infile is not None:
            raise MutantError("Infile reading not implemented for Carette mutant pool datasets!")
        # Just use the normal init, except self._mutants_by_position should contain Carette mutants, not normal ones
        Insertional_mutant_pool_dataset.__init__(self, cassette_end, reads_are_reverse, multi_dataset, infile)
        self._mutants_by_position = keybased_defaultdict(blank_Carette_mutant_with_pos)

    # TODO method for filling in genome-side Carette data!  Something based on add_alignment_reader_to_data too

    def add_Carette_genome_side_alignments_to_data(self, genome_side_infile, best_only=False, multiple_alignment=False, 
                                                   max_distance_for_best=MAX_POSITION_DISTANCE):
        # TODO docstring!
        if self.summary.cassette_end not in SEQ_ENDS:
            raise MutantError("Can't run add_Carette_genome_side_alignments_to_data if cassette_end is unknown!")
        if self.summary.reads_are_reverse not in [True,False]:
            raise MutantError("Can't run add_Carette_genome_side_alignments_to_data if reads_are_reverse is unknown!")
        # We'll be going over a lot of genome-side reads, which are the paired-end-sequencing mates of 
        #  the cassette-side reads used for mutant positions. 

        ### First make a dictionary with read IDs as keys and mutant positions as values, 
        # so we can quickly identify which read goes with which mutant 
        # (some/most reads won't go with any mutants, if their cassette-side reads were unaligned or discarded)
        read_ID_to_position = {}
        for mutant in self:
            for read_ID in mutant.read_IDs:
                read_ID_to_position[read_ID] = mutant.position

        ### Go over all the reads in the HTSeq_alignment_reader, add them to the appropriate mutant
        if genome_side_infile.endswith('.sam'):
            read_IDs = set()
            for aln in HTSeq.SAM_Reader(genome_side_infile):
                read_ID = aln.read.name
                # grab the mutant that has the matching read; if no mutant matches the read ID, go on to next read
                try:                mutant_pos = read_ID_to_position[read_ID]
                except KeyError:    continue
                mutant = self.get_mutant(mutant_pos)
                # for multiple alignments, just give the sequence and ignore the position entirely. 
                #  (have to keep track of read IDs, since there are likely to be multiple lines per read!)
                # TODO implement those properly, with giving the alignment positions?  TRICKY to figure out how to print them!
                if multiple_alignment:
                    if best_only:
                        raise MutantError("Multiple alignments make no sense with best_only flag!")
                    if read_ID in read_IDs:
                        continue
                    mutant.add_Carette_read(new_position='MULTIPLE', HTSeq_alignment=None, seq=aln.read.seq)
                # for normal alignments, use full alignment to add read info
                else:
                    if read_ID in read_IDs:
                        raise MutantError("File %s wasn't supposed to be multiple alignments - why is read %s there twice??"%(
                            genome_side_infile, read_ID))
                    # convert the genome-side read HTSeq alignment data to an "insertion" position of the end of it!
                    genome_side_position = get_Carette_pos_from_flanking_region_pos(aln.iv, self.summary.cassette_end, 
                                                                                    self.summary.reads_are_reverse)
                    if best_only:   mutant.improve_best_Carette_read(genome_side_position, aln, max_distance_for_best)
                    else:           mutant.add_Carette_read(genome_side_position, aln)
                read_IDs.add(read_ID)
        # For fasta files, just assume it's unaligned, or multiple if that's set, and use the sequences.
        elif genome_side_infile.endswith('.fa'):
            if best_only:
                raise MutantError("Fasta files (multiple or unaligned) make no sense with best_only flag!")
            if not multiple_alignment:
                print "Assuming reads in %s are unaligned."%genome_side_infile
                position_val = 'UNALIGNED'
            else:
                position_val = 'MULTIPLE'
            for name,seq in parse_fasta(genome_side_infile):
                try:                mutant_pos = read_ID_to_position[name]
                except KeyError:    continue
                mutant = self.get_mutant(mutant_pos)
                mutant.add_Carette_read(new_position=position_val, HTSeq_alignment=None, seq=seq)
        else:
            raise MutantError("Don't know how to deal with file %s - not a .fa or .sam file!"%genome_side_infile)

    def read_data_from_file(self, infile, assume_new_sequences=False):
        raise MutantError("Not implemented for Carette datasets!")

    def _set_merged_genome_info(self, gene_annotation_header_values, total_genes_in_genome_values):
        raise MutantError("Not implemented for Carette datasets!")

    def populate_multi_dataset(self, source_dataset_dict, overwrite=False, check_gene_data=True):
        raise MutantError("Not implemented for Carette datasets!")

    def merge_other_dataset(self, other_dataset):
        raise MutantError("Not implemented for Carette datasets!")
        # MAYBE-TODO this MIGHT actually work for Carette datasets, but I'm not sure and don't want to spend the time checking

    def find_genes_for_mutants(self, genefile, detailed_features=False, N_run_groups=3, verbosity_level=1):
        """ To each mutant in the dataset, add the gene it's in (look up gene positions for each mutant using genefile).

        ALSO add gene data to all the Carette-genome-side-read mutants inside each mutant!

        If detailed_features is True, also look up whether the mutant is in an exon/intron/UTR.
        Read the file in N_run_groups passes to avoid using up too much memory/CPU.
        """ 
        # Group all the mutants by chromosome, so that I can go over each chromosome in genefile separately
        #   instead of reading in all the data at once (which uses a lot of memory)
        #  Inclue both the main mutants, AND all the Carette genome-side read sub-mutants!
        mutants_by_chromosome = defaultdict(set)
        for mutant in self:
            mutants_by_chromosome[mutant.position.chromosome].add(mutant)
            for Carette_read_data in mutant.Carette_genome_side_reads:
                try:                    mutants_by_chromosome[Carette_read_data.position.chromosome].add(Carette_read_data)
                except AttributeError:  continue
        self._find_genes_for_mutant_list(mutants_by_chromosome, genefile, detailed_features, N_run_groups, verbosity_level)

    def add_gene_annotation(self, annotation_file, if_standard_Phytozome_file=None, custom_header=None, print_info=False):
        """ Add gene annotation to each mutant, based on annotation_file. See parse_gene_annotation_file doc for detail."""
        gene_annotation_dict = self._get_gene_annotation_dict(annotation_file, if_standard_Phytozome_file, custom_header, print_info)
        # add the annotation info to each mutant (or nothing, if gene has no annotation), 
        #  and to the Carette data sub-mutants
        N_annotated = 0
        for mutant in self:
            for curr_mutant in [mutant] + mutant.Carette_genome_side_reads:
                if_annotated = self._add_gene_annotation_to_mutant(curr_mutant, gene_annotation_dict)
                if if_annotated:    N_annotated += 1
        if print_info:          print "Added %s annotations"%N_annotated
        elif not N_annotated:   print "Warning: No gene annotations found!"
        # LATER-TODO add this to the gene-info run-test case!

    def print_detailed_Carette_data(self, OUTPUT=sys.stdout, sort_data_by=None, max_distance=MAX_POSITION_DISTANCE):
        """ Write detailed Carette data (all reads per mutant) to separate file.
        """
        # TODO docstring!

        # TODO should probably add header
        # TODO add annotation!

        ### sort all mutants by position or readcount (for single datasets only for now), or don't sort at all
        sorted_mutants = self._sort_data(sort_data_by)

        ### Quick summary
        N_total = len(self)
        # using sum because you can't do len on generators
        N_single_genomic = sum(1 for m in self if m.Carette_N_distinct_positions(max_distance)[0]==1)
        N_single_chrom = sum(1 for m in self if m.Carette_N_genomic_chromosomes==1)
        N_cassette = sum(1 for m in self if m.Carette_N_distinct_positions(max_distance)[1]>0)
        OUTPUT.write("# %s mutants total; %s have one genomic location; "%(N_total, 
                                                                           value_and_percentages(N_single_genomic, [N_total]))
                     +"%s have locations on only one genomic chromosome; "%(value_and_percentages(N_single_chrom, [N_total]))
                     +"%s also have one or more cassette locations.\n"%(value_and_percentages(N_cassette, [N_total])) )

        # TODO add info on how many actually have ANY genomic-side reads!  And out of those, how many are confirmed vs not.

        ### Print data for all mutants
        for mutant in sorted_mutants:
            mutant.Carette_print_detail(OUTPUT, max_distance)

    # TODO need to overwrite or modify print_summary (could probably use refactoring!)

    # TODO need to overwrite or modify print_data, to add columns that give some information on the Carette data!


############################### Further analysis/plotting utilities #####################################

# TODO should those be somewhere else?  In deconvolution, or plotting utilities, or separate Carette_utilities file?

def add_Carette_to_deconvolution_data(DECONV_header, DECONV_data, Carette_5prime, Carette_3prime, 
                                      max_allowed_distance, min_weird_distance):
    """ Add Carette data (two extra fields) to deconvolution output mutant list.
    """
    # TODO more detail in docstring!
    NEW_deconv_data = []
    NEW_deconv_header = DECONV_header[:2] + 'Carette_status Carette_distance'.split() + DECONV_header[2:]
    status_counts = defaultdict(int)
    for fields in DECONV_data:
        position = Insertion_position(fields[7], fields[8], full_position=fields[10], immutable=True)
        side = fields[6]
        if side=="5'":      mutant = Carette_5prime.get_mutant(position)
        else:               mutant = Carette_3prime.get_mutant(position)
        status = mutant.Carette_confirmation_status(side, max_allowed_distance, min_weird_distance=min_weird_distance)
        status_distance = mutant.Carette_max_confirmed_distance(max_allowed_distance)
        if isnan(status_distance): status_distance = '-'
        else:                            status_distance = str(status_distance)
        status_counts[status] += 1
        NEW_deconv_data.append(fields[:2] + tuple([status, status_distance]) + fields[2:])
    for status,count in sorted(status_counts.items()):
       print "%s: %s"%(status, value_and_percentages(count, [len(NEW_deconv_data)]))
    return NEW_deconv_header, NEW_deconv_data, status_counts
     # TODO unit-test!

def get_all_distances(deconv_data):
    """ Return lists of max confirmed distances for all mutants and separated into various categories.
    """
    # TODO more detail in docstring!
    distances_all, distances_genomic, distances_cassette = [], [], []
    distances_by_side_genomic, distances_by_quality_genomic = defaultdict(list), defaultdict(list)
    distances_by_side_cassette, distances_by_quality_cassette = defaultdict(list), defaultdict(list)
    for line in deconv_data:
        # ignore cases with no max confirmed distance
        if line[3] == '-':    continue
        dist = int(line[3])
        distances_all.append(dist)
        if is_cassette_chromosome(line[9]):
            distances_cassette.append(dist)
            distances_by_side, distances_by_quality = distances_by_side_cassette, distances_by_quality_cassette
        else:
            distances_genomic.append(dist)
            distances_by_side, distances_by_quality = distances_by_side_genomic, distances_by_quality_genomic
        distances_by_side[line[8]].append(dist)
        categories = line[4:6]
        if 'unmapped' in categories:    distances_by_quality['unmapped'].append(dist)
        elif 'poor' in categories:      distances_by_quality['poor'].append(dist)
        elif 'decent' in categories:    distances_by_quality['decent'].append(dist)
        elif 'good' in categories:      distances_by_quality['good'].append(dist)
        elif 'best' in categories:      distances_by_quality['best'].append(dist)
        else:                           print "weird categories??? %s"%' '.join(line)
    return distances_all, distances_genomic, distances_cassette, distances_by_side_genomic, distances_by_quality_genomic, distances_by_side_cassette, distances_by_quality_cassette
    # TODO unit-test!

def distance_histogram(distance_datasets, labels=None, colors=None, linestyles=None, N_bins=50, histtype='step', normalized=True, 
                       outfile=None, filetypes='svg png', x_max=None, sample_info=''):
    """ Plot histogram of the distances in the given datasets; optionally save to file. 
    """
    mplt.figure()
    max_y = 0
    for N, distances in enumerate(distance_datasets):
        if not len(distances):
            print "Dataset %s contains no data! Skipping."%N
            continue
        plot_kwargs = {}
        if labels:      plot_kwargs['label'] = "%s (%s)"%(labels[N], len(distances))
        if colors:      plot_kwargs['color'] = colors[N]
        if linestyles:  plot_kwargs['linestyle'] = linestyles[N]
        y_values, _, _ = mplt.hist(distances, bins=N_bins, histtype=histtype, normed=str(normalized), **plot_kwargs)
        max_y = max(max_y, max(y_values))
    mplt.xlabel('highest distance between cassette-side and genome-side reads')
    mplt.ylim(-max_y/20, max_y*1.1)
    if x_max is not None:
        mplt.xlim(mplt.xlim()[0], x_max)
    if normalized:  mplt.ylabel('fraction of mutants')
    else:           mplt.ylabel('number of mutants')
    if labels:
        mplt.legend()
    mplt.title("Carette confirmed distance distribution" + "\n(%s)"%sample_info if sample_info else "")
    if outfile:
        plotting_utilities.savefig(outfile, filetypes)
    mplt.close()

def plot_confirmed_dist_vs_percent(dataset, dataset_name, min_genomic_reads=1, min_conf_reads=1, max_allowed_distance=3000,
                                   markersize=4, alpha=0.3, color='black'):
    """ Make a max-confirmed-distance vs %-confirming-reads scatterplot.
    """
    D = max_allowed_distance
    basic_mutant_data = [(m.Carette_max_confirmed_distance(D), m.Carette_N_confirming_reads(D), m.Carette_N_non_confirming_reads(D)) 
                         for m in dataset]
    filtered_data = [x for x in basic_mutant_data if x[1]+x[2]>min_genomic_reads and x[1]>min_conf_reads]
    print "Filtering data (only include mutants with at least %s total genomic reads and at least %s confirming reads): %s/%s passed filter"%(min_genomic_reads, min_conf_reads, len(filtered_data), len(basic_mutant_data))
    max_conf_len = [x[0] for x in filtered_data]
    percent_wrong_reads = [x[2]/(x[1]+x[2])*100 for x in filtered_data]
    mplt.plot(max_conf_len, percent_wrong_reads, 
              marker='.', markerfacecolor=color, markersize=markersize, linestyle='None', markeredgecolor='None', alpha=alpha)
    mplt.xlabel('max distance between cassette-side sequence and matching genome-side read')
    mplt.ylabel("% reads that don't match the cassette-side sequence")
    mplt.ylim(-1, 101)
    mplt.title('Data on Carette genome-side reads confirming the cassette-side sequence,\n dataset %s, min %s total and %s confirming reads.'%(dataset_name, min_genomic_reads, min_conf_reads))


############################################### Tests ##############################################################

class Testing_all(unittest.TestCase):
    """ Unit-tests for position-related classes and functions. """

    def test__get_Carette_pos_from_flanking_region_pos(self):
        pos_tuple = ('C', 3, 7, '+')
        # note: For Carette cassette-side reads, 5' data has 5prime and reads_are_reverse=True, 3' data has 3prime and reads_are_reverse=False.  NOT DEALING WITH other cases - they should raise an error!
        self.assertRaises(MutantError, get_Carette_pos_from_flanking_region_pos, pos_tuple, '5prime', reads_are_reverse=False)
        self.assertRaises(MutantError, get_Carette_pos_from_flanking_region_pos, pos_tuple, '3prime', reads_are_reverse=True)
        ### case with +strand original read ([|||] is cassette, ---...--- is the paired-end read, * is desired position): 
        #   +strand cassette if 5' data, or -strand if 3':  *----......-----[|||||]
        pos_tuple = ('C', 3, 7, '+')
        pos = get_Carette_pos_from_flanking_region_pos(pos_tuple, '5prime', reads_are_reverse=True)
        assert (pos.chromosome == 'C' and pos.strand=='+' and pos.min_position==3)
        pos = get_Carette_pos_from_flanking_region_pos(pos_tuple, '3prime', reads_are_reverse=False)
        assert (pos.chromosome == 'C' and pos.strand=='-' and pos.min_position==3)
        ### the other case:  [|||||]----......-----*
        #   -strand cassette if 5' data, or +strand if 3'
        pos_tuple = ('C', 3, 7, '-')
        pos = get_Carette_pos_from_flanking_region_pos(pos_tuple, '5prime', reads_are_reverse=True)
        assert (pos.chromosome == 'C' and pos.strand=='-' and pos.min_position==7)
        pos = get_Carette_pos_from_flanking_region_pos(pos_tuple, '3prime', reads_are_reverse=False)
        assert (pos.chromosome == 'C' and pos.strand=='+' and pos.min_position==7)

    @staticmethod
    def _make_pos(pos):
        """ Convenience function to make the args to Fake_HTSeq_genomic_pos from an Insertion_position """
        return pos.chromosome, pos.strand, pos.min_position, pos.min_position+20

    def test__add_Carette_read(self):
        """ Also tests Carette_max_confirmed_distance """
        Fake_HTSeq_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment
        pos0 = Insertion_position('chr1','+',position_before=100)
        aln0 = Fake_HTSeq_alignment(seq='AAA', readname='read1', pos=self._make_pos(pos0))
        mutant = Insertional_mutant_Carette(pos0)
        assert mutant.total_read_count == 0
        assert isnan(mutant.Carette_max_confirmed_distance(10))
        mutant.add_read(aln0)
        assert mutant.total_read_count == 1
        pos1 = Insertion_position('chr1','+',position_before=0)
        aln1 = Fake_HTSeq_alignment(seq='AAA', readname='read1', pos=self._make_pos(pos1))
        mutant.add_Carette_read(pos1, HTSeq_alignment=aln1)
        pos2 = Insertion_position('chr2','+',position_before=0)
        aln2 = Fake_HTSeq_alignment(seq='AGA', readname='read1', pos=self._make_pos(pos2))
        mutant.add_Carette_read(pos2, HTSeq_alignment=aln2)
        pos3 = Insertion_position('chr1','-',position_before=0)
        aln3 = Fake_HTSeq_alignment(seq='ACA', readname='read1', pos=self._make_pos(pos3))
        mutant.add_Carette_read(pos3, HTSeq_alignment=aln3)
        assert len(mutant.Carette_genome_side_reads) == 3
        assert [m.position.chromosome for m in mutant.Carette_genome_side_reads] == 'chr1 chr2 chr1'.split()
        assert [m.position.strand for m in mutant.Carette_genome_side_reads] == '+ + -'.split()
        assert mutant.Carette_max_confirmed_distance(1000) == 100
        assert mutant.Carette_max_confirmed_distance(100) == 100
        assert isnan(mutant.Carette_max_confirmed_distance(10))

    def test__improve_best_Carette_read(self):
        """ Also tests Carette_max_confirmed_distance """
        Fake_HTSeq_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment
        # make the cassette-side read
        pos0 = Insertion_position('chr1','+',position_before=100)
        aln0 = Fake_HTSeq_alignment(seq='AAA', readname='read1', pos=self._make_pos(pos0))
        mutant = Insertional_mutant_Carette(pos0)
        assert mutant.position.min_position == 100
        assert mutant.total_read_count == 0
        assert isnan(mutant.Carette_max_confirmed_distance(1000))
        mutant.add_read(aln0)
        assert mutant.total_read_count == 1
        assert isnan(mutant.Carette_max_confirmed_distance(1000))
        # add bad read first, then better ones, make sure they're replaced (or not, if both are equally bad)
        pos1 = Insertion_position('chr2','+',position_before=10)
        aln1 = Fake_HTSeq_alignment(seq='AGA', readname='read1', pos=self._make_pos(pos1))
        mutant.improve_best_Carette_read(pos1, HTSeq_alignment=aln1, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 10
        assert isnan(mutant.Carette_max_confirmed_distance(1000))
        pos2 = Insertion_position('chr1','-',position_before=11)
        aln2 = Fake_HTSeq_alignment(seq='ACA', readname='read1', pos=self._make_pos(pos2))
        mutant.improve_best_Carette_read(pos2, HTSeq_alignment=aln2, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 10
        assert isnan(mutant.Carette_max_confirmed_distance(1000))
        pos3 = Insertion_position('chr1','+',position_before=12)
        aln3 = Fake_HTSeq_alignment(seq='AAA', readname='read1', pos=self._make_pos(pos3))
        mutant.improve_best_Carette_read(pos3, HTSeq_alignment=aln3, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 10
        assert isnan(mutant.Carette_max_confirmed_distance(1000))
        pos4 = Insertion_position('chr1','+',position_before=80)
        aln4 = Fake_HTSeq_alignment(seq='AAA', readname='read1', pos=self._make_pos(pos4))
        mutant.improve_best_Carette_read(pos4, HTSeq_alignment=aln4, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 80
        assert mutant.Carette_max_confirmed_distance(1000) == 20
        pos5 = Insertion_position('chr1','+',position_before=70)
        aln5 = Fake_HTSeq_alignment(seq='AAA', readname='read1', pos=self._make_pos(pos5))
        mutant.improve_best_Carette_read(pos5, HTSeq_alignment=aln5, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 70
        assert mutant.Carette_max_confirmed_distance(1000) == 30
        # add good read first, then worse ones, make sure it's NOT replaced (remember to start with new mutant!)
        mutant = Insertional_mutant_Carette(pos0)
        mutant.add_read(aln0)
        mutant.improve_best_Carette_read(pos5, HTSeq_alignment=aln5, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 70
        assert mutant.Carette_max_confirmed_distance(1000) == 30
        mutant.improve_best_Carette_read(pos1, HTSeq_alignment=aln1, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 70
        assert mutant.Carette_max_confirmed_distance(1000) == 30
        mutant.improve_best_Carette_read(pos2, HTSeq_alignment=aln2, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 70
        assert mutant.Carette_max_confirmed_distance(1000) == 30
        mutant.improve_best_Carette_read(pos3, HTSeq_alignment=aln3, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 70
        assert mutant.Carette_max_confirmed_distance(1000) == 30
        mutant.improve_best_Carette_read(pos4, HTSeq_alignment=aln4, max_distance=50)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 70
        assert mutant.Carette_max_confirmed_distance(1000) == 30
        # make sure wrong-strand closer-distance new read won't replace old right-strand one! (this was a bug once)
        pos0 = Insertion_position('chr1','+',position_before=100)
        aln0 = Fake_HTSeq_alignment(seq='AAA', readname='read1', pos=self._make_pos(pos0))
        mutant = Insertional_mutant_Carette(pos0)
        mutant.add_read(aln0)
        pos5 = Insertion_position('chr1','+',position_before=70)
        aln5 = Fake_HTSeq_alignment(seq='AAA', readname='read1', pos=self._make_pos(pos5))
        mutant.improve_best_Carette_read(pos5, HTSeq_alignment=aln5, max_distance=1000)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 70
        assert mutant.Carette_max_confirmed_distance(1000) == 30
        pos6 = Insertion_position('chr1','-',position_before=11)
        aln6 = Fake_HTSeq_alignment(seq='ACA', readname='read1', pos=self._make_pos(pos6))
        mutant.improve_best_Carette_read(pos6, HTSeq_alignment=aln6, max_distance=1000)
        assert len(mutant.Carette_genome_side_reads) == 1
        assert mutant.Carette_genome_side_reads[0].position.min_position == 70
        assert mutant.Carette_max_confirmed_distance(1000) == 30

if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    # if I wanted more control I could do this instead:
    #import os
    #unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromName(os.path.splitext(sys.argv[0])[0]))
    #   (autodetection of all tests - see http://docs.python.org/library/unittest.html#unittest.TestLoader)
    # there's probably also some way to easily get all tests from the current file without passing the name, but I haven't found it yet...
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main(argv=[sys.argv[0]])
