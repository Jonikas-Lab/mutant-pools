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

######################################## Calculating mappability over genome #######################################

### MEMORY USAGE ISSUES
# Keeping the full set of mappable/unmappable genome positions, in whatever form, takes a lot of memory!
# Use sys.getsizeof to get the sizes of basic objects:  
#  - an int or bool takes 24 bytes; 
#  - a 21-string takes 58, a 1-string takes 38
#  - sets and dicts take 1-2x as much memory as the objects themselves.
#  - lists actually take a lot less!  About 9 bits per int or char.
# The genome's about 113M bases, and we already know that most of it's mappable, i.e. unique sequences.  
#  - that many ints in a set (say 1.5x the cost) will be 4G!  
#  - storing that many 21-strings in a set (say 1.5x the cost) will be 9G!  
#  - plus even more for the seq:pos dict, depending on details of pos format...

### MAYBE-TODO Improvement options for memory issues:
# For the second pass and returned value: Instead of a list of genome positions, keep a genome-sized list filled with 0s and 1s?  With 0/1 as normal ints that would still be about 1G, but I'm sure switching to a more efficient representation (bitstring.BitArray?) would make it a LOT smaller. 
# For the first pass - NOT SURE!  How can I more efficiently store sets of strings?  MAYBE-TODO think about it more - it's currently the part that takes the most memory. 

## Issues with python not releasing memory, or something
# 
# Technically, on the whole genome, genome_mappable_slices should take 5-8GB, but RELEASE it once done (see "del" calls) - the output is only <1G.  But python seems to keep the memory, even if the output is about the size I expected!  So if I do two runs, say 20bp and 21bp flanks, I end up exceeding my 16G and going into swap... 
# genome_mappable_insertion_sites* doesn't take as much memory to run, though still some; the output is actually 2x bigger than genome_mappable_slices for genome_mappable_insertion_sites, since there are two mappable insertion positions per slice, and 2Nx bigger for genome_mappable_insertion_sites_multi (4x for two slice lengths, etc - comes out to about 3G). 
#
# So it's probably a good idea to pickle the output (after running genome_mappable_slices or genome_mappable_insertion_sites*), close python to release the memory (even though it'll take a LONG TIME, especially if there was swap space used), and open a new instance to unpickle the data and work on it.  
# BUT if I pickle the output of any of those, and unpickle it in a new interactive python shell, again the output itself doesn't end up too huge (see sizes above), BUT the unpickling seems to take more memory and not release it again! Ugh. So maybe I really should convert it to bitarrays or something. 


def genome_mappable_slices(slice_len, genome_seq=None, print_info=True):
    """ Return a chrom:position_list dict giving the positions of all the unique slice_len seq slices. 
    
    Look at all slices of length slice_len in each chromosome (overlapping - for AATGG with len 2, look at AA, AT, TG, GG)
     return the positions of only the ones that showed up exactly once (in forward or reverse-complement - 
      sequences that have a reverse-complement version later in the genome, or that are palindromes, don't show up).
    Positions are the 1-based position of the first base of each slice (i.e. in AATTGGCC, the position of AAT is 1).

    Genome_seq can be a chr_name:chr_seq dict, or the path to a fasta file containing the genome; default file will be used if None. 

    CAVEAT: if you run this on the whole Chlamy genome, around 100Mb, it can take 5-8G of RAM, and returns a 1G data structure!
    """
    # keep original genome_seq, because we'll want it for the second pass, and generators only work once
    original_genome_seq = genome_seq
    if genome_seq is None:           genome_seq = os.path.expanduser(DEFAULT_ALL_GENOME_FILE)
    if isinstance(genome_seq, str):  genome_seq = basic_seq_utilities.parse_fasta(genome_seq)
    else:                            genome_seq = genome_seq.iteritems()
    ### This is done in two passes over the whole genome, to improve memory usage. The original version was done in one pass.
    ### First pass - go over all chromosome and get a set of unique and non-unique sequences. 
    sequences_seen_once = set()
    sequences_seen_twice = set()
    chrom_names = set()
    N_total_slices = 0
    if print_info:  print "PASS 1 - finding unique sequences (printing only chromosomes over 1Mb):"
    for chrom_name,chrom_seq in genome_seq:
        if print_info and len(chrom_seq)>1000000: 
            print "  %s (%s)...  %s"%(chrom_name, basic_seq_utilities.format_base_distance(len(chrom_seq)), time.ctime())
        if chrom_name in chrom_names:
            raise Exception("%s seen twice in genome_seq!"%(chrom_name))
        chrom_names.add(chrom_name)
        # go over each slice_len-sized sequence slice and check if it's unique
        for slice_start,slice_seq in basic_seq_utilities.generate_seq_slices(chrom_seq, slice_len, step=1):
            N_total_slices += 1
            # look at both the forward and the reverse-complement sequence of the current slice, 
            #  but only store the forward (+strand) versions (no point in storing both if we're RC-ing each current sequence)
            slice_seq_RC = basic_seq_utilities.reverse_complement(slice_seq)
            # if a sequence is its own reverse-complement, it cannot be uniquely mappable - ignore
            #  (no point in even adding it to any dictionaries, since any future copies of that sequence will be caught here too)
            if slice_seq == slice_seq_RC:
                continue
            # if this is the first time slice_seq shows up, save it, and save its position
            #  (note - originally I did "{slice_seq,RC} & sequences_seen_once | sequences_seen_twice" here, HORRIBLY SLOW)
            if not (slice_seq in sequences_seen_once or slice_seq_RC in sequences_seen_once 
                    or slice_seq in sequences_seen_twice or slice_seq_RC in sequences_seen_twice):
                sequences_seen_once.add(slice_seq)
            # if this is the second time slice_seq shows up, note it as non-unique, 
            #  and remove both it AND its RC from the unique set and the unique position dict.
            #   (use set.discard(x) to avoid raising an error when one of them isn't there)
            elif (slice_seq in sequences_seen_once or slice_seq_RC in sequences_seen_once):
                sequences_seen_twice.add(slice_seq)
                sequences_seen_once.discard(slice_seq)
                sequences_seen_once.discard(slice_seq_RC)
            # if slice_seq is already known to be non-unique (in seen-twice but not seen-once), ignore
    # this is a HUGE data structure, release it as soon as possible - MAYBE-TODO not sure if that works...
    del sequences_seen_twice    

    ### Second pass - go over all chromosomes again, and save the positions of known unique sequences.
    # restart genome_seq generator
    genome_seq = original_genome_seq
    if genome_seq is None:              genome_seq = os.path.expanduser(DEFAULT_ALL_GENOME_FILE)
    if isinstance(genome_seq, str):     genome_seq = basic_seq_utilities.parse_fasta(genome_seq)
    else:                               genome_seq = genome_seq.iteritems()
    unique_seq_positions_by_chrom = {}
    if print_info:  print "PASS 2 - getting uniquely mappable positions (printing only chromosomes over 1Mb):"
    for chrom_name,chrom_seq in genome_seq:
        if print_info and len(chrom_seq)>1000000: 
            print "  %s (%s)...  %s"%(chrom_name, basic_seq_utilities.format_base_distance(len(chrom_seq)), time.ctime())
        unique_seq_positions_by_chrom[chrom_name] = []
        # go over each slice_len-sized sequence slice and check if it's unique
        for slice_start,slice_seq in basic_seq_utilities.generate_seq_slices(chrom_seq, slice_len, step=1):
            # we already have the full set of unique sequences - now just store the position in unique_seq_positions_by_chrom
            #  if it's a unique sequence.
            if slice_seq in sequences_seen_once:
                unique_seq_positions_by_chrom[chrom_name].append(slice_start)
    # this is a HUGE data structure, release it as soon as possible - MAYBE-TODO not sure if that works...
    del sequences_seen_once

    if print_info:
        N_mappable_slices = sum([len(pos_list) for pos_list in unique_seq_positions_by_chrom.values()])
        fraction_mappable = N_mappable_slices/N_total_slices if N_total_slices else float('NaN')
        print("%.0f%% of %sbp slices are mappable (%s out of %s total on %s chromosomes)"
              %(fraction_mappable*100, slice_len, N_mappable_slices, N_total_slices, len(chrom_names)))
    return unique_seq_positions_by_chrom


def genome_mappable_insertion_sites(flanking_region_length=21, mappable_slice_pos_dict=None, genome_seq=None, 
                                    include_strand=True, end_sequenced='5prime', print_info=True):
    """ Return all uniquely mappable genomic insertion sites, as dict with chrom or (chrom,strand) keys and pos_list values.

    The mappability is based on the mappability of the flanking region, using flanking_region_length, 
     on either side of the insertion (both sides are counted). 

    The actual positions are the positions of the base BEFORE the insertion; the position lists are sorted. 
    For example if the position 100-120 flanking region is mappable, that means an insertion at 99-100 is mappable, 
     and one at 120-121 is mappable, depending on the orientation: 
      if the flanking region is 5prime, then a +strand 120-121 and a -strand 99-100 insertion is mappable; 
      if 3prime, the strands are inverted. 
    Each dict value is a list of the base positions before the mappable insertion site (1-based).
    Each dict value (position_list) is sorted.  Each cassette orientation gets a separate entry in it.
    If include_strand is True, all positions are separated by strand, so the keys in the returned dictionary are
     (chromosome,strand) tuples - so if both a +strand and a -strand insertion at position X of chr1 is mappable, 
      position X will show up once in the ('chr1','+') list and once in ('chr1','-').
     If include_strand is False, the positions are merged, and the return dict keys are just chromosome names,
      so in the same example, position X will just show up twice on the 'chr1' list; 
      if position X is mappable only in one orientation, it'll show up once in the list.
     (The False option is meant for plotting, e.g. with mutant_plotting_utilities.mutant_positions_and_data)

    If mappable_slice_pos_dict is not None, it'll be assumed to be the output of genome_mappable_slices (with slice_len 
     the same as flanking_region_length), and genome_seq will be ignored; 
     if it's None, genome_mappable_slices will be used to generate that data (takes a while, and a lot of memory!).

    Genome_seq can be a list of (chr_name,chr_seq) tuples, or the path to a fasta file containing the genome; 
     default file will be used if None. 
    """
    if end_sequenced not in basic_seq_utilities.SEQ_ENDS:
        raise Exception('Invalid end_sequenced %s! Must be in %s.'%(end_sequenced, basic_seq_utilities.SEQ_ENDS))
    reads_are_reverse = False if end_sequenced=='5prime' else True
    # get a chromosome:mappable_slice_start_pos_list dict, if not provided already (this is the only use of genome_seq)
    if mappable_slice_pos_dict is None:
        mappable_slice_pos_dict = genome_mappable_slices(flanking_region_length, genome_seq, print_info)
    # convert the mappable slice chromosome/position values into mappable insertion locations, 
    #  treating them as end_sequenced flanking regions - each flanking region should give TWO mappable positions, 
    #   one on each side, with opposite strands (with the insertion position strand depending on end_sequenced)
    mappable_position_data = defaultdict(list)
    for chrom,pos_list in mappable_slice_pos_dict.iteritems():
        if print_info and len(pos_list)>500000: 
            print "  %s (%s mappable slices)...  %s"%(chrom, len(pos_list), time.ctime())
        for pos in pos_list:
            for strand in basic_seq_utilities.SEQ_STRANDS:
                flanking_region_pos = (chrom, pos, pos+flanking_region_length-1, strand)
                # the easiest way of getting the right insertion position is to just use the same function I normally use
                #  for making actual Insertion_position objects from sequenced flanking region position data
                insertion_pos_object = mutant_analysis_classes.get_insertion_pos_from_flanking_region_pos(flanking_region_pos, 
                                                                                              end_sequenced, reads_are_reverse)
                if include_strand:
                    mappable_position_data[chrom,insertion_pos_object.strand].append((insertion_pos_object.min_position))
                else:
                    mappable_position_data[chrom].append(insertion_pos_object.min_position)
    return general_utilities.sort_lists_inside_dict(mappable_position_data)


def genome_mappable_insertion_sites_multi(flanking_region_lengths=[20,21], mappable_slice_pos_dicts=None, genome_seq=None, 
                                          include_strand=True, end_sequenced='5prime', print_info=True):
    """ Run genome_mappable_insertion_sites multiple times, with different flank length and optionally mappable dict values.

    Flanking_region_lengths should be a list of flanking_region_length values; 
     mappable_slice_pos_dicts can be None, or a same-length list of values. 
    For each pair of those values, genome_mappable_insertion_sites will be run with the matching arguments. 
    The final output will include the sum of position lists from all runs. 

    See genome_mappable_insertion_sites docstring for more info.
    """
    if mappable_slice_pos_dicts == None:    mappable_slice_pos_dicts = [None for _ in flanking_region_lengths]
    full_mappable_position_data = defaultdict(list)
    # for each set of inputs, run genome_mappable_insertion_sites to generate new data, 
    #  and then merge the new data into full_mappable_position_data
    for flanking_region_length,mappable_slice_pos_dict in zip(flanking_region_lengths, mappable_slice_pos_dicts):
        # get data for 
        curr_mappable_position_data = genome_mappable_insertion_sites(flanking_region_length, mappable_slice_pos_dict, genome_seq, include_strand, end_sequenced, print_info)
        for chrom,new_pos_list in curr_mappable_position_data.iteritems():
            full_mappable_position_data[chrom].extend(new_pos_list)
    # sort the final positions (also converts from defaultdict to dict)
    return general_utilities.sort_lists_inside_dict(full_mappable_position_data)


################################# Simulating random size-N dataset, various options #######################################

def get_20bp_fraction(dataset):
    """ Return the fraction of mutants in dataset that only have 20bp flanking regions (no 21bp ones).

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list of mutants (mutant_analysis_classes.Insertional_mutant instances).
    """
    max_lengths = [max([len(s) for s in mutant.sequences_and_counts.keys()]) for mutant in dataset]
    assert max_lengths.count(20) + max_lengths.count(21) == len(dataset), "Some sequences are outside the 20-21bp range!"
    return max_lengths.count(20)/len(dataset)
    # LATER-TODO should this be here, or a mutant_analysis_classes.Insertional_mutant_pool_dataset method or something?


### simulate dataset with N randomly positioned mutants, taking into account mappable/unmappable positions!
# TODO actual simulation function!  Use get_20bp_fraction to decide how many of the simulated mutants should be 20bp vs 21bp. 


### simulate dataset with N randomly positioned mutants, matching the gap-size distribution of the real dataset, and taking into account mappable/unmappable positions

# LATER-TODO
# MAYBE-TODO would it be possible to also match the real hotspots and coldspots of the real dataset?


################################# Randomly chosen subsets of real dataset #######################################

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
        #  - arguments to genome_mappable_insertion_sites are (flanking_region_length=21, mappable_slice_pos_dict=None, 
        #           genome_seq=None, include_strand=True, end_sequenced='5prime', print_info=True):

        ### only testing whether there are 0 unique sequences, or a non-zero number
        def _test_all_empty(slice_len, genome):
            """ Test whether all the genome_mappable_* functions return empty sets/lists/dicts/whatever. """
            def _test_empty(D): return sum([len(l) for l in D.values()])==0
            outcomes = set()
            # test genome_mappable_slices
            slice_data = genome_mappable_slices(slice_len, genome_seq=genome, print_info=False)
            outcomes.add(_test_empty(slice_data))
            # test genome_mappable_insertion_sites and genome_mappable_insertion_sites_multi in all variants
            for include_strand in (True,False):
                for end in basic_seq_utilities.SEQ_ENDS:
                    args = (include_strand,end, False)
                    outcomes.add(_test_empty(genome_mappable_insertion_sites(slice_len, slice_data, None, *args)))
                    outcomes.add(_test_empty(genome_mappable_insertion_sites(slice_len, None, genome, *args)))
                    for N in (1,2,5):
                        outcomes.add(_test_empty(genome_mappable_insertion_sites_multi([slice_len]*N, [slice_data]*N, None, *args)))
                        outcomes.add(_test_empty(genome_mappable_insertion_sites_multi([slice_len]*N, None, genome, *args)))
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
        # help functions for testing all the functions in parallel
        def _read_raw_data(raw_data):
            """ get a dict of lists of ints or (int.str) tuples from a raw string. """
            formatted_data = defaultdict(list)
            for chr_data in raw_data.split(', '):
                chrom, pos_data = chr_data.strip().split(': ')
                # pos_data can be ints (1 25 301) or ints with strand info (1- 25+ 301+)
                for x in pos_data.split(' '):
                    if x[-1] in '+-':   formatted_data[chrom,x[-1]].append(int(x[:-1]))
                    else:               formatted_data[chrom].append(int(x))
            for val in formatted_data.values():
                val.sort()
            return dict(formatted_data)
        def _test_all(slice_len, genome, raw_slice_data, raw_pos_data_5prime):
            """ Test genome_mappable_slices and all variants of genome_mappable_insertion_sites* against expected output. """
            # get the expected output data from the simplified string formats
            slice_data = _read_raw_data(raw_slice_data)
            pos_data_5prime = _read_raw_data(raw_pos_data_5prime)
            # from pos_data_5prime, make pos_data_3prime (just switch all strands) and pos_data_no_strand (remove strand info)
            pos_data_3prime, pos_data_no_strand = {}, defaultdict(list)
            for (chrom,strand),pos_list in pos_data_5prime.items():
                pos_data_3prime[chrom, '+' if strand=='-' else '-'] = pos_list
                pos_data_no_strand[chrom] += pos_list
                pos_data_no_strand[chrom].sort()
            # check genome_mappable_slices output (and save it for later)
            new_slice_data = genome_mappable_slices(slice_len, genome, False)
            assert new_slice_data == slice_data
            # now try running genome_mappable_insertion_sites with both the raw slice_len/genome data, and the new_slice_data; 
            for include_strand,end,pos_data in ((True,'5prime',pos_data_5prime), (True,'3prime',pos_data_3prime), 
                                                (False,'5prime',pos_data_no_strand), (False,'3prime',pos_data_no_strand)):
                    args = (include_strand,end, False)
                    assert genome_mappable_insertion_sites(slice_len, slice_data, None, *args) == pos_data
                    assert genome_mappable_insertion_sites(slice_len, None, genome, *args) == pos_data
                    assert genome_mappable_insertion_sites_multi([slice_len], None, genome, *args) == pos_data
                    for extra in ([], [{}], [{}, {}, {}]):
                        assert genome_mappable_insertion_sites_multi([slice_len]*(len(extra)+1), [slice_data]+extra, None, 
                                                                     *args) == pos_data
        # how genome_mappable_insertion_sites works if end is 5prime: 
        #   a mappable 1-2 flanking region means a +strand insertion at 2-? and a -strand insertion at ?-1 will be mappable;
        #  if end is 3prime, it's the same positions but opposite strands.

        # a single-base-repeat chromosome with slice size equal to chromosome size has all unique sequences;
        #  same for two such chromosomes that aren't reverse-complement
        _test_all(3, {'a':'AAA'},                  'a: 1',        'a: 0- 3+')
        _test_all(3, {'a':'AAA', 'b':'GGG'},       'a: 1, b: 1',  'a: 0- 3+, b: 0- 3+')
        # same test, but with the genome read from a fasta file
        _test_all(3, 'test_data/INPUT_genome0.fa', 'a: 1, b: 1',  'a: 0- 3+, b: 0- 3+')

        # a single chromosome that's not a palindrome or a base-repeat has all unique sequences
        #  (or just some if slice_len is 1 and some bases show up twice)
        # (note - the first case here has two different-strand mappable insertions in one position, 
        #  so the genome_mappable_insertion_sites has one position repeated twice - important to check that!)
        _test_all(1, {'a':'AG'},   'a: 1 2',  'a: 0- 1+ 1- 2+') 
        _test_all(2, {'a':'AAT'},  'a: 1',    'a: 0- 2+')
        curr_genome = {'a':'AGA'}
        _test_all(1, curr_genome,  'a: 2',    'a: 1- 2+')
        _test_all(2, curr_genome,  'a: 1 2',  'a: 0- 2+ 1- 3+')

        # test genome_mappable_insertion_sites_multi on the two curr_genome cases added together, different flank sizes
        #  (getting the data directly from genome, or from two genome_mappable_slices results)
        #  (only testing the no-strand version for simplicity - MAYBE-TODO test all versions?)
        assert genome_mappable_insertion_sites_multi([1,2], None, curr_genome, False, '5prime', False) == {'a':[0,1,1,2,2,3]}
        slices_1 = genome_mappable_slices(1, curr_genome, False)
        slices_2 = genome_mappable_slices(2, curr_genome, False)
        for fl_both,slices_both in (([1,2],[slices_1,slices_2]), ([2,1],[slices_2,slices_1])):
            assert genome_mappable_insertion_sites_multi(fl_both, slices_both, None, False, '5prime', False) == {'a':[0,1,1,2,2,3]}

        
    # LATER-TODO add unit-tests for other stuff!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
