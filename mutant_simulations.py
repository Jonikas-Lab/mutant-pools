#! /usr/bin/env python
"""
Plotting utilities specifically for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import unittest
import time
import os
import math
import random
from collections import defaultdict
# other packages
import numpy
import scipy
# my modules
import mutant_analysis_classes
import general_utilities
import basic_seq_utilities

DEFAULT_NUCLEAR_GENOME_FILE = '~/experiments/reference_data/genomes_and_indexes/Chlre4-nm.fa'
DEFAULT_ALL_GENOME_FILE = '~/experiments/reference_data/genomes_and_indexes/Chlre4-nm_chl-mit.fa'
DEFAULT_GENOME_CASSETTE_FILE = '~/experiments/reference_data/genomes_and_indexes/Chlre4-nm_chl-mit_cassette-pMJ013b.fa'

########################################### Mutant simulation functions ##########################################

### mappability over genome

### MEMORY USAGE ISSUES
# Keeping the full set of mappable/unmappable genome positions, in whatever form, takes a lot of memory!
# Use sys.getsizeof to get the sizes of basic objects:  
#  - an int or bool takes 24 bytes; 
#  - a 21-string takes 58, a 1-string takes 38
#  - sets and dicts take 1-2x as much memory as the objects themselves.
#  - lists actually take a lot less!  About 9 bits per int or char.
# The genome's about 113M bases, and we already know that most of it's mappable, i.e. unique sequences.  
#  - storing that many 21-strings in a set (say 1.5x the cost) will be 9G!  
#  - that many ints in a set (say 1.5x the cost) will be 4G!  
#  - plus even more for the seq:pos dict, depending on details of pos format...

### improvement options:
# 
# Make unique_seq_positions values contain only the start position, not start and end, since they're all the same length
# Make unique_seq_positions values use a number for the chromosome instead of a name
# -> these to get the dictionary down to 100M * (58+24+24) * 1.5 - about 15G.  Still too much!
#
# Change to keeping a genome-sized list filled with 0s and 1s?  That would be 1G, which is doable...  EXCEPT that I can't actually do that - I need to know which sequence corresponds to which position, otherwise if later I find a second copy of a given sequence, I won't know which position to mark as non-unique for the first sequence!
#
# I could do two passes over the genome... First just get a set of unique sequences (which is <9G, hopefully), THEN go over all the sequences again, check which ones are unique, and full the genome-sized 0/1 list based on that.
#
# Note for the genome-sized 0/1 list - I actually want just 0/1, not full ints, so I'm sure switching to a more efficient representation (bitstring.BitArray?) would make it a LOT smaller. 

def genome_mappable_slices(slice_len, genome_seq=None, print_info=True):
    """ Return a sequence:position dict giving the positions of all the unique seq slices. 
    
    Look at all the slices of length slice_len in each chromosome, both forward and reverse-complement;
     return the positions of only the ones that showed up exactly once.
    Positions are (chrom, slice_start, slice_end) tuples; start/end positions are 1-based inclusive
     (i.e. in AATTGGCC, the position of AA is 1-2).
    Positions are always given on the +strand; obviously if a +strand sequence is unique, 
     the same start-end -strand sequence, which is its reverse-complement, is also unique. 
     However, both strands are checked - palindromic sequences won't show up as unique.

    Genome_seq can be a chr_name:chr_seq dict, or the path to a fasta file containing the genome; default file will be used if None. 
    """
    # TODO this takes a LOT of memory!!! What to do?
    if genome_seq is None:           genome_seq = os.path.expanduser(DEFAULT_ALL_GENOME_FILE)
    if isinstance(genome_seq, str):  genome_seq = basic_seq_utilities.parse_fasta(genome_seq)
    else:                            genome_seq = genome_seq.iteritems()
    sequences_seen_twice = set()
    unique_seq_positions = {}
    chrom_names = set()
    N_total_slices = 0
    # go over all chromosomes in genome file
    for chrom_name,chrom_seq in genome_seq:
        if print_info and len(chrom_seq)>1000000: 
            print "processing %s (%sbp)...  %s"%(chrom_name, len(chrom_seq), time.ctime())
        if chrom_name in chrom_names:
            raise Exception("%s seen twice in genome_seq!"%(chrom_name))
        chrom_names.add(chrom_name)
        # go over each slice_len-sized sequence slice and check if it's unique
        for slice_start,slice_seq in basic_seq_utilities.generate_seq_slices(chrom_seq, slice_len, step=1):
            N_total_slices += 1
            # look at both the forward and the reverse-complement sequence of the current slice, 
            #  but only store the forward (+strand) versions (no point in doing both)
            slice_seq_RC = basic_seq_utilities.reverse_complement(slice_seq)
            # if a sequence is its own reverse-complement, it cannot be uniquely mappable - ignore
            #  (no point in even adding it to any dictionaries, since any future copies of that sequence will be caught here too)
            if slice_seq == slice_seq_RC:
                continue
            # if this is the first time slice_seq shows up, save it, and save its position
            #  (note - originally I did "{slice_seq,RC} & unique_seq_positions | sequences_seen_twice" here, HORRIBLY SLOW)
            if not (slice_seq in unique_seq_positions or slice_seq_RC in unique_seq_positions 
                    or slice_seq in sequences_seen_twice or slice_seq_RC in sequences_seen_twice):
                unique_seq_positions[slice_seq] = (chrom_name, slice_start, slice_start+slice_len-1)
            # if this is the second time slice_seq shows up, note it as non-unique, 
            #  and remove both it AND its RC from the unique set and the unique position dict.
            #   (use set.discard(x) and dict.pop(x,None) to avoid raising an error when one of them isn't there)
            elif (slice_seq in unique_seq_positions or slice_seq_RC in unique_seq_positions):
                sequences_seen_twice.add(slice_seq)
                unique_seq_positions.pop(slice_seq,None)
                unique_seq_positions.pop(slice_seq_RC,None)
            # if slice_seq is already known to be non-unique (in seen-twice but not seen-once), ignore
    if print_info:
        N_mappable_slices = len(unique_seq_positions)
        fraction_mappable = N_mappable_slices/N_total_slices if N_total_slices else float('NaN')
        print("%s%% of %sbp slices are mappable (%s out of %s total, counting forward and reverse, on %s chromosomes)"
              %(fraction_mappable*100, slice_len, N_mappable_slices, N_total_slices, len(chrom_names)))
    return unique_seq_positions


def genome_mappable_insertion_positions(flanking_region_length=21, genome_seq=None, end_sequenced='5prime', print_info=True):
    """ Return a set of all uniquely mappable genomic insertion position OBJECTS, given flanking_region_length. 

    Insertion sites are mutant_analysis_classes.Insertion_position instances (with only one end known);
     end_sequenced must be '5prime' or '3prime', specifying which end of the insertion has the mapped flanking region.
    Give the length of the flanking regions - sometimes 21bp ones will be mappable but 20bp ones won't!

    Genome_seq can be a list of (chr_name,chr_seq) tuples, or the path to a fasta file containing the genome; 
     default file will be used if None. 
    """
    # get a list of unique seq positions for all mappable slices
    unique_flanking_region_positions = genome_mappable_slices(flanking_region_length, genome_seq, print_info).values()
    # convert those into mutant_analysis_classes.Insertion_position objects - 
    #  there will be two insertion_positions for each slice, one on each end.
    ### convert the flanking region (chrom,start_pos,end_pos) list to a list of Insertion_position instances:
    #   each flanking region should give TWO mappable positions, one on each side, with opposite strands
    #    (with the insertion position strand depending on end_sequenced)
    unique_insertion_positions = set()
    if end_sequenced not in basic_seq_utilities.SEQ_ENDS:
        raise Exception('Invalid end_sequenced %s! Must be in %s.'%(end_sequenced, basic_seq_utilities.SEQ_ENDS))
    reads_are_reverse = False if end_sequenced=='5prime' else True
    for flanking_region_pos in unique_flanking_region_positions:
        for strand in basic_seq_utilities.SEQ_STRANDS:
            unique_insertion_positions.add(mutant_analysis_classes.get_insertion_pos_from_flanking_region_pos(
                                                            flanking_region_pos + (strand,), end_sequenced, reads_are_reverse))
    return unique_insertion_positions


def genome_mappable_insertion_sites_by_chrom(flanking_region_lengths=[20,21], genome_seq=None, print_info=True,
                                             position_sets=None):
    """ Return chromosome:position_list dict with all uniquely mappable genomic insertion SITES, omitting strand, for plotting.

    Based on genome_mappable_insertion_positions (see docstring for that!), but transformed to a chromosome:position_list dict, 
     where the positions are the min_position of each Insertion_position instance.  Each position_list is sorted.
    Strand information is discarded; if a given insertion position is mappable regardless of insertion orientation, 
     it'll show up on the list twice; if it's mappable only in one orientation, it'll show up once.
    This is to make it usable as other_dataset argument element for mutant_plotting_utilities.mutant_positions_and_data.

    Two ways of providing input:
     - provide flanking_region_lengths and genome_seq - genome_mappable_insertion_positions will be run
        for each value in flanking_region_lengths
     - provide position_sets - a list of outputs from genome_mappable_insertion_positions to be used directly 
        to avoid re-running; flanking_region_lengths and genome_seq are ignored

    If a single value is given as flanking_region_lengths, just run genome_mappable_insertion_positions with that; 
     if a list of values is given, run genome_mappable_insertion_positions on each and add the results together.

    There's no end_sequenced arg like in genome_mappable_insertion_positions - 5prime is always used,
     since only strand information changes based on end_sequenced, and we're ignoring strand information anyway.
    """
    chromosome_position_lists = defaultdict(list)
    # if position_sets are given (must be sets of mutant_analysis_classes.Insertion_position instances), 
    #  just grab the min_position value of each Insertion_position, adding to appropriate chromosome list.
    if position_sets is not None:
        for position_set in position_sets:
            for mappable_pos in position_set:
                chromosome_position_lists[mappable_pos.chromosome].append(mappable_pos.min_position)
    # if position_sets not given, generate them based on flanking_region_lengths and genome_seq.
    else:
        if isinstance(flanking_region_lengths, int):
            flanking_region_lengths = [flanking_region_lengths]
        for fl_len in flanking_region_lengths:
            for mappable_pos in genome_mappable_insertion_positions(fl_len, genome_seq, '5prime', print_info):
                chromosome_position_lists[mappable_pos.chromosome].append(mappable_pos.min_position)
    # sort the positions!
    return general_utilities.sort_lists_inside_dict(chromosome_position_lists)


### simulate dataset with N randomly positioned mutants, taking into account mappable/unmappable positions!

def get_20bp_fraction(dataset):
    """ Return the fraction of mutants in dataset that only have 20bp flanking regions (no 21bp ones).

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list of mutants (mutant_analysis_classes.Insertional_mutant instances).
    """
    max_lengths = [max([len(s) for s in mutant.sequences_and_counts.keys()]) for mutant in dataset]
    assert max_lengths.count(20) + max_lengths.count(21) == len(dataset), "Some sequences are outside the 20-21bp range!"
    return max_lengths.count(20)/len(dataset)
    # LATER-TODO should this be here, or a mutant_analysis_classes.Insertional_mutant_pool_dataset method or something?


# TODO actual simulation function!  Use get_20bp_fraction to decide how many of the simulated mutants should be 20bp vs 21bp. 


### simulate dataset with N randomly positioned mutants, matching the gap-size distribution of the real dataset, and taking into account mappable/unmappable positions

# LATER-TODO
# MAYBE-TODO would it be possible to also match the real hotspots and coldspots of the real dataset?


### number of genes with 1+/2+/etc mutants vs number of mutants (randomly chosen mutant subsets)

def _genes_with_N_mutants(mutants, N):
    """ Given a mutant list, return the number of genes with at least N mutants on the list. """
    gene_mutant_counts = defaultdict(int)
    for m in mutants:
        if m.gene not in mutant_analysis_classes.SPECIAL_GENE_CODES.all_codes:
            gene_mutant_counts[m.gene] += 1
    return len([1 for count in gene_mutant_counts.values() if count>=N])
   
def _genes_with_no_mutants(mutants, total_genes=17114):
    """ Given a list of mutants, return the number of genes with no mutants on the list. 
    The total number of genes can be given as an argument; default is 17114 (from Phytozome chlamy v4.3).
    """
    return total_genes - genes_with_N_mutants(mutants, 1)

   
def gene_counts_for_mutant_subsets(dataset, max_N_mutants=3, step_size=100, single_subset_size=None):
    """ Return numbers of genes with N mutants for different-sized random subsets of dataset, or a single subset size.

    Return an N:list_of_gene_numbers dictionary, where N is each value between 1 and max_N_mutants, 
     and list_of_gene_numbers contains the number of genes with at least N mutants 
      in randomly chosen subsets of dataset, of sizes starting at 0 and going up to the full dataset size in step_size steps.
     (So if dataset contains 200 mutants, and step_size is 100, each list will have 3 values, for 0, 100 and 200 mutants.)
     (Note that the last value is not for the full dataset size, but for the closest lower number divisible by step_size.)
    If single_subset_size is not None, ignore step_size, and just return an N:gene_number dictionary 
     for a randomly chosen subset of size single_subset_size.

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list of mutants (mutant_analysis_classes.Insertional_mutant instances).
    """
    # extract list of mutants from dataset
    if isinstance(dataset, mutant_analysis_classes.Insertional_mutant_pool_dataset):    mutants = list(dataset)
    else:                                                                               mutants = dataset
    # get subset_sizes list based on either single_subset_size or step_size
    if single_subset_size is not None:
        subset_sizes = [single_subset_size]
        step_size = 0
    else:
        subset_sizes = range(0, len(mutants), step_size)
    if step_size > len(dataset) or single_subset_size > len(dataset):
        raise ValueError("step_size and single_subset_size can't be greater than dataset size!")
    # only shuffle the dataset if we're not going to be using the full one anyway
    if not step_size == len(dataset) or single_subset_size == len(dataset):
        random.shuffle(mutants)
    # get the gene counts for each mutant number
    gene_counts = {}
    for N_mutants in range(1,max_N_mutants+1):
        gene_counts[N_mutants] = []
        for subset_size in subset_sizes:
            gene_counts[N_mutants].append(_genes_with_N_mutants(mutants[:subset_size], N_mutants))
    return gene_counts


################################################# Testing etc ##################################################

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__genome_mappable_everything(self):
        # testing three related functions in parallel, since they use the same cases
        #  - arguments to genome_mappable_slicesare (slice_len, genome_seq, print_info=True)
        #  - arguments to genome_mappable_insertion_positions are 
        #       (flanking_region_length=21, genome_seq=None, end_sequenced='5prime', print_info=True)
        #  - arguments to genome_mappable_insertion_sites_by_chrom are (flanking_region_length=21, genome_seq=None, print_info=True)

        ### only testing whether there are 0 unique sequences, or a non-zero number
        def _test_all_empty(slice_len, genome):
            """ Test whether all the genome_mappable_* functions return empty sets/lists/dicts/whatever. """
            outcomes = set([genome_mappable_slices(slice_len, genome, False) == {}])
            for end in basic_seq_utilities.SEQ_ENDS:
                outcomes.add(genome_mappable_insertion_positions(slice_len, genome, end, False) == set())
            outcomes.add(genome_mappable_insertion_sites_by_chrom(slice_len, genome, False) == {})
            if len(outcomes) != 1:
                raise Exception("Inconsistent outcomes in _test_all_empty in test__genome_mappable_everything!")
            return outcomes.pop()
        # no unique sequences regardless of slice_len: 
        #  empty genome, two identical or reverse-complement chromosomes, one palindromic chromosome, 
        for slice_len in (1,2,3,4):
            assert _test_all_empty(slice_len, {})
            assert _test_all_empty(slice_len, {'a':'AAA', 'b':'AAA'})
            assert _test_all_empty(slice_len, {'a':'AAA', 'b':'TTT'})
            assert _test_all_empty(slice_len, {'a':'AATT'})
        # one chromosome with 3+ tandem repeats - no unique sequences as long as slice_len <= repeat_len
        #  (this doesn't apply for 2 tandem repeats, because in ACTACT, CTA is still unique! But in ACTACTACT, it's present twice.)
        for repeat_seq in ('ATC', 'AATCCG', 'ACATGACGAGACGGG'):
            for N_repeats in (3,4,5,10):
                chrom_seq = repeat_seq * N_repeats
                for slice_len in range(1, len(repeat_seq)+2):
                    assert _test_all_empty(slice_len, {'a':chrom_seq})
        # no unique sequences if slice_len<chromosome_len:
        #  one chromosome that has repeating substrings or all matching substrings with own reverse-complement
        for slice_len in (1,2):
            assert _test_all_empty(slice_len, {'a':'AAA'})
            assert _test_all_empty(slice_len, {'a':'ATA'})
        # no unique sequences ONLY IF slice_len==1:  any case with more than one AT or GC base.
        for slice_len,if_empty in ((1,True), (2,False), (3,False), (10,False)):
            assert _test_all_empty(slice_len, {'a':'AA'}) == if_empty
            assert _test_all_empty(slice_len, {'a':'ATT'}) == if_empty
            assert _test_all_empty(slice_len, {'a':'GCC'}) == if_empty
            assert _test_all_empty(slice_len, {'a':'G', 'b':'GC'}) == if_empty
            assert _test_all_empty(slice_len, {'a':'AG', 'b':'TC'}) == if_empty
        # cases that should have non-zero unique sequences
        assert not _test_all_empty(3, {'a':'AAA'})
        assert not _test_all_empty(3, {'a':'AAA', 'b':'GGG'})
        assert not _test_all_empty(1, {'a':'AGA'})
        assert not _test_all_empty(2, {'a':'AGA'})
        assert not _test_all_empty(2, {'a':'AAT'})

        ### more detailed testing of actual non-empty outputs

        # help functions for testing genome_mappable_insertion_positions, which has a complicated output format
        def _I(*args, **kwargs):
            return mutant_analysis_classes.Insertion_position(*args, immutable=True, **kwargs)
        def _test_P(slice_len, genome, pos_list_5prime):
            pos_list_5prime = [tuple(single.split(' ')) for single in pos_list_5prime.split(',  ')]
            assert genome_mappable_insertion_positions(slice_len, genome, '5prime', False) == set([_I(*a) for a in pos_list_5prime])
            pos_list_3prime = [(c, '+' if s=='-' else '-', p) for (c,s,p) in pos_list_5prime]
            assert genome_mappable_insertion_positions(slice_len, genome, '3prime', False) == set([_I(*a) for a in pos_list_3prime])
        # how genome_mappable_insertion_positions works if end is 5prime: 
        #   a mappable 1-2 flanking region means a +strand insertion at 2-? and a -strand insertion at ?-1 will be mappable;
        #  if end is 3prime, it's the same positions but opposite strands.

        # a single-base-repeat chromosome with slice size equal to chromosome size has all unique sequences;
        #  same for two such chromosomes that aren't reverse-complement
        assert genome_mappable_slices(3, {'a':'AAA'}, False) == {'AAA':('a',1,3)}
        _test_P(3, {'a':'AAA'}, 'a - ?-1,  a + 3-?')
        assert genome_mappable_insertion_sites_by_chrom(3, {'a':'AAA'}, False) == {'a':[0,3]}

        assert genome_mappable_slices(3, {'a':'AAA', 'b':'GGG'}, False) == {'AAA':('a',1,3), 'GGG':('b',1,3)}
        _test_P(3, {'a':'AAA',  'b':'GGG'}, 'a - ?-1,  a + 3-?,  b - ?-1,  b + 3-?')
        assert genome_mappable_insertion_sites_by_chrom(3, {'a':'AAA', 'b':'GGG'}, False) == {ch:[0,3] for ch in 'ab'}

        # a single chromosome that's not a palindrome or a base-repeat has all unique sequences
        #  (or just some if slice_len is 1 and some bases show up twice)
        # (note - the first case here has two different-strand mappable insertions in one position, 
        #  so the genome_mappable_insertion_sites_by_chrom has one position repeated twice - important to check that!)

        assert genome_mappable_slices(1, {'a':'AG'}, False) == {'A':('a',1,1), 'G':('a',2,2)}
        _test_P(1, {'a':'AG'}, 'a - ?-1,  a + 1-?,  a - ?-2,  a + 2-?')
        assert genome_mappable_insertion_sites_by_chrom(1, {'a':'AG'}, False) == {'a':[0,1,1,2]}

        assert genome_mappable_slices(1, {'a':'AGA'}, False) == {'G':('a',2,2)}
        _test_P(1, {'a':'AGA'}, 'a - ?-2,  a + 2-?')
        assert genome_mappable_insertion_sites_by_chrom(1, {'a':'AGA'}, False) == {'a':[1,2]}
        assert genome_mappable_slices(2, {'a':'AGA'}, False) == {'AG':('a',1,2), 'GA':('a',2,3)}
        _test_P(2, {'a':'AGA'}, 'a - ?-1,  a + 2-?,  a - ?-2,  a + 3-?')
        assert genome_mappable_insertion_sites_by_chrom(2, {'a':'AGA'}, False) == {'a':[0,1,2,3]}
        # same genome but with both flank lengths (results are added together)
        assert genome_mappable_insertion_sites_by_chrom([1,2], {'a':'AGA'}, False) == {'a':[0,1,1,2,2,3]}
        # also test that genome_mappable_insertion_sites_by_chrom works properly when provided position_sets directly:
        pos_set_slice1 = genome_mappable_insertion_positions(1, {'a':'AGA'}, '5prime', False)
        pos_set_slice2 = genome_mappable_insertion_positions(2, {'a':'AGA'}, '5prime', False)
        assert genome_mappable_insertion_sites_by_chrom(position_sets=[pos_set_slice1], print_info=False) == {'a':[1,2]}
        assert genome_mappable_insertion_sites_by_chrom(position_sets=[pos_set_slice2], print_info=False) == {'a':[0,1,2,3]}
        for pos_slice_both in ([pos_set_slice1,pos_set_slice2], [pos_set_slice2,pos_set_slice1]):
            assert genome_mappable_insertion_sites_by_chrom(position_sets=pos_slice_both, print_info=False) == {'a':[0,1,1,2,2,3]}

        assert genome_mappable_slices(2, {'a':'AAT'}, False) == {'AA':('a',1,2)}
        _test_P(2, {'a':'AAT'}, 'a - ?-1,  a + 2-?')
        assert genome_mappable_insertion_sites_by_chrom(2, {'a':'AAT'}, False) == {'a':[0,2]}

        
    # LATER-TODO add unit-tests for other stuff!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
