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
from mutant_analysis_classes import MutantError, Insertion_position, Insertional_mutant, Insertional_mutant_pool_dataset
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

    def __init__(self, *args, **kwargs):
        """ Same as for Insertional_mutant, plus add empty set of read IDs. """
        Insertional_mutant.__init__(self, *args, **kwargs)
        self.read_IDs = set()
        self.Carette_other_side_reads = []

    def add_read(self, HTSeq_alignment, read_count=1, treat_unknown_as_match=False, dataset_name=None):
        """ Add a read to the data (or multiple identical reads, if read_count>1); return True if perfect match.

        Runs Insertional_mutant.add_read, plus adds read ID to self.read_IDs.
        """
        Insertional_mutant.add_read(self, HTSeq_alignment, read_count, treat_unknown_as_match, dataset_name)
        self.read_IDs.add(HTSeq_alignment.read.name)

    def merge_mutant(self, *args, **kwargs):
        raise MutantError("merge_mutant NOT IMPLEMENTED on multi-dataset mutant object!")
        # TODO implement?  Just need to merge the read IDs and Carette read data...

    def add_Carette_read(self, HTSeq_alignment, read_count=1, treat_unknown_as_match=False):
        pass
        # TODO implement!

    # TODO add new method that gives some kind of simple summary of the Carette data - all same chromosome, different, etc
    # TODO add new method that infers the approximate real insertion location!

    # TODO TODO TODO finish implementing!  What other methods to add/overwrite?

# TODO what about mutants that combine normal MmeI and Carette data - are those going to be any different?


class Insertional_mutant_pool_dataset_Carette(Insertional_mutant_pool_dataset):
    # TODO docstring!

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

    # TODO overwrite add_alignment_reader_to_data to include Carette data?  Or can that just be done for mutants?
    # TODO method for filling in genome-side Carette data!  Something based on add_alignment_reader_to_data too

    def add_alignment_reader_to_data(self, HTSeq_alignment_reader, uncollapse_read_counts=False, 
                                     ignore_cassette=False, cassette_only=False, treat_unknown_as_match=False):
        """ Adds all alignments from the reader to the mutant data; currently based only on position, but that may change. 

        Input must be a list/generator/etc of HTSeq.Alignment objects (usually an HTSeq.SAM_Reader).

        Set uncollapse_read_counts to True if the original deepseq data was collapsed to unique sequences using
         fastx_uncollapser before alignment, to get the correct original read counts.

        Treat_unknown_as_match governs whether alignments with no detailed information are treated as perfect or not.

        Different ways of treating cassette reads:
         - by default (if all *cassette* args are false) treat cassette reads normally
         - if ignore_cassette==True, ignore cassette reads in the data and list them as removed in the header
         - if cassette_only==True, ignore all OTHER reads and only include cassette reads!
        """
        # LATER-TODO actually instead of cassette_only it might be good to just generate two separate mutant-sets, normal and cassette, with an option called separate_cassette or something, and print them to separate files - but that's more complicated, and right now I don't have the setup for a single dataset having multiple mutant-sets (although I guess I will have to eventually, for removed mutants etc). Right now I do it in mutant_count_alignments.py, which works but there's a lot of code repetition...
        if self.multi_dataset:  raise MutantError("add_alignment_reader_to_data not implemented for multi-datasets!")
        if ignore_cassette and cassette_only:
            raise MutantError("Only one of ignore_cassette and cassette_only arguments can be True - mutually exclusive!")

        summ = self.summary
        if summ.cassette_end == '?':
            raise MutantError("Cannot add data from an alignment reader if cassette_end isn't specified! Please set the "
          +"summary.cassette_end attribute of this Insertional_mutant_pool_dataset instance to one of %s first."%SEQ_ENDS)
        if summ.reads_are_reverse == '?':
            raise MutantError("Cannot add data from an alignment reader if reads_are_reverse isn't set! Please set the "
                  +"reads_are_reverse attribute of this Insertional_mutant_pool_dataset instance to True/False first.")

        # Go over all the reads in the HTSeq_alignment_reader, add them to dataset
        for aln in HTSeq_alignment_reader:
            if uncollapse_read_counts:      read_count = get_seq_count_from_collapsed_header(aln.read.name)
            else:                           read_count = 1
            # if read is unaligned, add to unaligned count and skip to the next read
            if (not aln.aligned) or (aln.iv is None):
                summ.non_aligned_read_count += read_count
                continue
            # get the cassette insertion position (as an Insertion_position object) - USE IMMUTABLE POSITIONS BY DEFAULT
            position = get_insertion_pos_from_flanking_region_pos(aln.iv, summ.cassette_end, summ.reads_are_reverse, 
                                                                  immutable_position=True)
            # if read is aligned to cassette and should be ignored, add to the right count and skip to the next read
            if ignore_cassette and is_cassette_chromosome(position.chromosome):
                summ.ignored_region_read_counts[position.chromosome] += read_count
                continue
            # or if we only want cassette reads, skip non-cassette ones!
            elif cassette_only and not is_cassette_chromosome(position.chromosome):
                summ.ignored_region_read_counts['NON-CASSETTE'] += read_count
                continue
            # grab the right mutant based on the position, and add the reads to it; 
            curr_mutant = self.get_mutant(position)
            curr_mutant.add_read(aln, read_count, treat_unknown_as_match=treat_unknown_as_match)
        # special case for when we don't know the specific unaligned categories, but we know total non-aligned is 0, 
        #  so the specific categories must be 0 too:
        if summ.non_aligned_read_count==0:  summ.unaligned, summ.multiple_aligned = 0, 0

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

    # TODO need to overwrite or modify print_summary (could probably use refactoring!)

    # TODO need to overwrite or modify this, to add columns that give some information on the Carette data!
    def print_data(self, OUTPUT=sys.stdout, sort_data_by=None, N_sequences=None, header_line=True, header_prefix="# "):
        """ Print full data, one line per mutant: position data, gene info, read counts, optionally sequences.
        (see the file header line for exactly what all the output fields are).

        For a normal dataset, print position/gene info (7 fields), total_reads, perfect_reads, N_sequence_variants, 
         the N_sequences most common sequences and counts (alternating, taking 2*N_sequences fields), 
         and gene annotation if it exists (arbitrary number of tab-separated fields).

        For a multi-dataset, print position/gene info (7 fields), the main sequence (taken over all the datasets),
         then reads_in_<x> and perfect_in_<x> for each dataset (total and perfect reads - total 2*N_datasets fields), 
         then gene annotation if it exists.

        Data is printed to OUTPUT, which should be an open file object (stdout by default).
         Output is tab-separated, with optional header starting with "# ".  
        """
        if self.multi_dataset and N_sequences is not None and N_sequences!=1:
            raise MutantError("Only one sequence can currently be printed in print_data for multi-datasets!")
        if N_sequences is None and not self.multi_dataset:     N_sequences = 2

        ### print the header line (different for normal and multi-datasets)
        # MAYBE-TODO should printing the gene info be optional?  Maybe... 
        #   Would save space when there isn't any meaningful gene info to print anyway.
        if header_line:
            header = ['chromosome','strand','min_position','full_position', 'gene','orientation','feature']
            if not self.multi_dataset:
                header += ['total_reads','perfect_reads', 'N_sequence_variants']
                for N in range(1,N_sequences+1):    
                    header += ['read_sequence_%s'%N, 'seq_%s_count'%N]
            else:
                header += ['main_sequence']
                for dataset_name in self.dataset_order:
                    header += ['reads_in_%s'%dataset_name, 'perfect_in_%s'%dataset_name]
            header += self.gene_annotation_header
            OUTPUT.write(header_prefix + '\t'.join(header) + "\n")

        ### sort all mutants by position or readcount (for single datasets only for now), or don't sort at all
        all_mutants = iter(self)
        if sort_data_by=='position':
            data = sorted(all_mutants, key = lambda m: m.position)
            # x.position here is an Insertion_position object and has a sensible cmp function
        elif sort_data_by=='read_count':
            if self.multi_dataset:  
                raise MutantError("Sorting by readcount in print_data not implemented for multi-datasets!")
            data = sorted(all_mutants, key = lambda m: (m.total_read_count,m.perfect_read_count,m.position), reverse=True)
        else:
            data = all_mutants

        # create "empty" annotation line with the correct number of fields, for genes that weren't in annotation file
        if self.gene_annotation_header:
            missing_gene_annotation_data = ['NO GENE DATA'] + ['' for x in self.gene_annotation_header[:-1]]

        ### for each mutant, print the mutant data line (different for normal and multi-datasets)
        for mutant in data:
            mutant_data = [mutant.position.chromosome, mutant.position.strand, mutant.position.min_position, 
                           mutant.position.full_position, mutant.gene, mutant.orientation, mutant.gene_feature] 
            if not self.multi_dataset:
                mutant_data += [mutant.total_read_count, mutant.perfect_read_count, mutant.unique_sequence_count]
                for N in range(1,N_sequences+1):
                    mutant_data += list(mutant.get_main_sequence(N))
                    # MAYBE-TODO also give the length and number of mutations for each sequence? Optionally?  
                    #   Length is easy, but do I even keep track of mutation number?  I probably should...
            else:
                mutant_data += [mutant.get_main_sequence(1)[0]]
                for dataset_name in self.dataset_order:
                    mutant_data += [mutant.by_dataset[dataset_name].total_read_count, 
                                    mutant.by_dataset[dataset_name].perfect_read_count]
            # add gene annotation, or a line with the right number of fields if gene annotation is missing
            #  (or if the mutant has no such attribute at all - possible for older-format datasets)
            if self.gene_annotation_header:
                try:
                    if mutant.gene_annotation:  mutant_data += mutant.gene_annotation
                    else:                       mutant_data += missing_gene_annotation_data
                except AttributeError:          mutant_data += missing_gene_annotation_data
            OUTPUT.write('\t'.join([str(x) for x in mutant_data]))
            OUTPUT.write('\n')


############################################### Unit-tests ##############################################################

class Testing_all(unittest.TestCase):
    """ Unit-tests for position-related classes and functions. """
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
