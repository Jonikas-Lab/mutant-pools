#!/usr/bin/env python
"""
Module containing classes and functions for analysis of deepseq data related to insertional mutant libraries.

This is a module to be imported and used by other programs.  Running it directly runs the built-in unit-test suite.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2013
"""

from __future__ import division
# basic libraries
import sys, re
import unittest
from collections import defaultdict
# other libraries
import HTSeq
# my modules
from mutant_analysis_classes import MutantError, Insertion_position, Insertional_mutant, Insertional_mutant_pool_dataset, is_cassette_chromosome
from general_utilities import keybased_defaultdict, value_and_percentages
# there's a "from parse_annotation_file import parse_gene_annotation_file" in one function, not always needed


# TODO need to implement new Insertion_position variant with uncertainty!  And modify the gene/feature finding function to deal with that...


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

    def add_Carette_read(self, HTSeq_alignment, read_count=1):
        self.Carette_genome_side_reads.append(HTSeq_alignment)
        # TODO do we want to parse that further now, or just leave it like this until later?
        # TODO should probably get gene info for those!  And annotation info!

    def Carette_N_distinct_positions(self, max_distance=1000):
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
        # add all the genome-side read positions
        for read in self.Carette_genome_side_reads:
            if read.aligned:
                positions_by_chrom_strand[(read.iv.chrom, read.iv.strand)].append(read.iv.start)
        # count total number of dictinct positions - different chromosomes or strands, or distance > max_distance
        total_distinct_positions_genome = 0
        total_distinct_positions_cassette = 0
        # for each chromosome, go over all positions and only count ones every 1000bp as distinct
        for chrom_strand, positions in positions_by_chrom_strand.items():
            positions.sort()
            distinct_positions = [positions[0]]
            for pos in positions[1:]:
                if (pos-distinct_positions[-1]) > max_distance:
                    distinct_positions.append(pos)
            if is_cassette_chromosome(chrom_strand[0]):     total_distinct_positions_cassette += len(distinct_positions)
            else:                                           total_distinct_positions_genome += len(distinct_positions)
        N_genomic_chromosomes = len(set([c for c in zip(*positions_by_chrom_strand)[0] if not is_cassette_chromosome(c)]))
        return total_distinct_positions_genome, total_distinct_positions_cassette, N_genomic_chromosomes
    
    # TODO add new method that infers the approximate real insertion location!

    # TODO add new method that prints info!
    def Carette_print_detail(self, OUTPUT=sys.stdout, max_distance=1000):
        # TODO make max_distance default value a class variable, so it's the same in different functions?
        ### print summary line
        N_distinct_genome, N_distinct_cassette, N_chroms_genome = self.Carette_N_distinct_positions(max_distance)
        N_genomic_reads = sum(1 for aln in self.Carette_genome_side_reads if not is_cassette_chromosome(aln.iv.chrom))
        OUTPUT.write(" * %s genome-side reads for %s distinct genomic positions (on %s chromosomes), "%(N_genomic_reads, 
                                                                                        N_distinct_genome, N_chroms_genome)
                     +"plus %s distinct cassette positions:\n"%(N_distinct_cassette))
        ### print cassette-side position
        OUTPUT.write("Cassette-side position::: %s, %s reads, top seq %s (%s reads)\n"%((self.position, self.total_read_count) + 
                     self.get_main_sequence()))
        ### print lines for each of the Carette genome-side reads
        # sort the Carette reads by alignment position (chrom,strand,pos; unaligned go last)
        OUTPUT.write("Carette protocol genome-side reads::: (SAM format)\n")
        def _position_sort_key(HTSeq_read):
            if not HTSeq_read.aligned:
                return 1
            return (0, HTSeq_read.iv.chrom, HTSeq_read.iv.strand, HTSeq_read.iv.start)
        for alignment in sorted(self.Carette_genome_side_reads, key = _position_sort_key):
            OUTPUT.write(alignment.get_sam_line() + '\n')
            # TODO SAM format is stupid, switch to something sensible here! 
            #   make sure direction is same as from cassette-side read if they actually match
            #   add gene/annotation data for each position!

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

    def add_Carette_genome_side_alignments_to_data(self, HTSeq_alignment_reader):
        # TODO docstring!
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
        for aln in HTSeq_alignment_reader:
            # if read is unaligned, add to unaligned count and skip to the next read
            read_ID = aln.read.name
            # TODO do the two paired-ends have slightly different IDs?
            # grab the mutant that has the matching read; if no mutant matches the read ID, go on to next read
            try:                mutant_pos = read_ID_to_position[read_ID]
            except KeyError:    continue
            mutant = self.get_mutant(mutant_pos)
            mutant.add_Carette_read(aln)

    def read_data_from_file(self, infile, assume_new_sequences=False):
        raise MutantError("Not implemented for Carette datasets!")

    def _set_merged_genome_info(self, gene_annotation_header_values, total_genes_in_genome_values):
        raise MutantError("Not implemented for Carette datasets!")

    def populate_multi_dataset(self, source_dataset_dict, overwrite=False, check_gene_data=True):
        raise MutantError("Not implemented for Carette datasets!")

    def merge_other_dataset(self, other_dataset):
        raise MutantError("Not implemented for Carette datasets!")
        # MAYBE-TODO this MIGHT actually work for Carette datasets, but I'm not sure and don't want to spend the time checking

    # TODO need to write new function to write detailed Carette data (all Carette reads per mutant) to separate file!
    def print_detailed_Carette_data(self, OUTPUT=sys.stdout, sort_data_by=None, max_distance=1000):
        # TODO docstring!

        ### sort all mutants by position or readcount (for single datasets only for now), or don't sort at all
        sorted_mutants = self._sort_data(sort_data_by)

        ### Quick summary
        N_total = len(self)
        # using sum because you can't do len on generators
        N_single_genomic = sum(1 for m in self if m.Carette_N_distinct_positions(max_distance)[0]==1)
        N_single_chrom = sum(1 for m in self if m.Carette_N_distinct_positions(max_distance)[2]==1)
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


############################################### Unit-tests ##############################################################

class Testing_all(unittest.TestCase):
    """ Unit-tests for position-related classes and functions. """

    def test__(self):
        print "NO UNIT-TESTS IMPLEMENTED YET!"


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    # if I wanted more control I could do this instead:
    #import os
    #unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromName(os.path.splitext(sys.argv[0])[0]))
    #   (autodetection of all tests - see http://docs.python.org/library/unittest.html#unittest.TestLoader)
    # there's probably also some way to easily get all tests from the current file without passing the name, but I haven't found it yet...
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main(argv=[sys.argv[0]])
