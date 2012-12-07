#! /usr/bin/env python
"""
Simulation functions for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import unittest
import time
import os
import math
import random
import collections
# other packages
import numpy
import scipy
import scipy.stats
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
# my modules
import general_utilities
import basic_seq_utilities
import mutant_analysis_classes
from mutant_utilities import get_histogram_data_from_positions, get_mutant_positions_from_dataset, get_chromosome_lengths, get_20bp_fraction

######################################## Calculating mappability over genome #######################################

### MEMORY USAGE ISSUES
    # Keeping the full set of mappable/unmappable genome positions, in whatever form, takes a lot of memory!
    # Use sys.getsizeof to get the sizes of basic objects:  
    #  - an int or bool takes 24 bytes; 
    #  - a 21-string takes 58, a 1-string takes 38
    #  - sets and dicts take 1-2x as much memory as the objects themselves.
    #  - lists actually take a lot less!  About 9 bits per int or char (TODO is that really true? Geoffrey says not. Did I test it with big ints and non-identical chars?).
    #  - numpy arrays take a lot less, but I think sys.getsizeof lies about them because they're objects... Still, in practice on real mappability data they seem to take about 2x less than lists, which gets it down to something manageable.  They're also about 2x bigger than the non-numpy version when pickled (1.3G vs 800M), but that's not as big an issue, I'd say.
    # 
    # The genome's about 113M bases, and we already know that most of it's mappable, i.e. unique sequences.  
    #  - that many ints in a set (say 1.5x the cost) will be 4G!  
    #  - storing that many 21-strings in a set (say 1.5x the cost) will be 9G!  
    #  - plus even more for the seq:pos dict, depending on details of pos format...

    ### MAYBE-TODO Improvement options for memory issues:
    # 
    # For the second pass and returned value: Instead of a list of genome positions, keep a genome-sized list filled with 0s and 1s?  With 0/1 as normal ints that would still be about 1G, but I'm sure switching to a more efficient representation (bitstring.BitArray?) would make it a LOT smaller.  Actually it might be smaller as a numpy array too, those may account for element size...
    # 
    # Could I do everything as numpy arrays instead of making lists (with append/extend) and then converting to numpy arrays?  NOT SURE.  I could write the list-creation as a generator, but you can't make numpy arrays from generators without converting to a list first, I checked...
    # 
    # For the first pass - NOT SURE!  How can I more efficiently store sets of strings?  One thing I could do would be use ints or bitarrays or something instead of strings, since all my strings are actually just made of ACTG, with 2 bits per character, not the full range of characters (well, there may be Ns, but I can skip those...). 
    # MAYBE-TODO think about it more - it's currently the part that takes the most memory.  But I only run it once, so that's not really a big issue - it's mostly the returned values that are important.

    ## Issues with python not releasing memory, or something
    # 
    # Technically, on the whole genome, genome_mappable_slices should take 5-8GB, but RELEASE it once done (see "del" calls) - the output is only <1G.  But python seems to keep the memory, even if the output is about the size I expected!  So if I do two runs, say 20bp and 21bp flanks, I end up exceeding my 16G and going into swap... 
    # genome_mappable_insertion_sites* doesn't take as much memory to run, though still some; the output is actually 2x bigger than genome_mappable_slices for genome_mappable_insertion_sites, since there are two mappable insertion positions per slice, and 2Nx bigger for genome_mappable_insertion_sites_multi (4x for two slice lengths, etc - comes out to about 3G). 
    #
    # So it's probably a good idea to pickle the output (after running genome_mappable_slices or genome_mappable_insertion_sites*), close python to release the memory (even though it'll take a LONG TIME, especially if there was swap space used), and open a new instance to unpickle the data and work on it.  
    # BUT if I pickle the output of any of those, and unpickle it in a new interactive python shell, again the output itself doesn't end up too huge (see sizes above), BUT the unpickling seems to take more memory and not release it again! Ugh. So maybe I really should convert it to bitarrays or something.  Although numpy arrays helped some. 


def genome_mappable_slices(slice_len, genome_seq=None, print_info=True):
    """ Return a chrom:position_array dict giving the positions of all the unique slice_len seq slices. 
    
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
            # make sure everything's uppercase
            slice_seq = slice_seq.upper()
            # if there are Ns in the sequence, just consider it unmappable
            #  MAYBE-TODO do something more complicated?  Technically sequences with one N might be unique,
            #  line ANA if there are no AAA, AGA, ACA or ATA in the whole genome...  Ignoring that for now.  
            #  I'm not even sure if bowtie ever considers genome sequences with Ns mappable...
            if 'N' in slice_seq:
                continue
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

    # Convert the result to a numpy array, it's way faster!  But I have to convert it afterward rather than doing it
    #  as a numpy array to start with, because you can't append to numpy arrays.
    unique_seq_positions_by_chrom = {chrom: numpy.array(pos_list) for (chrom, pos_list) in unique_seq_positions_by_chrom.items()}
    return unique_seq_positions_by_chrom


def genome_mappable_insertion_sites(flanking_region_length=21, mappable_slice_pos_dict=None, genome_seq=None, 
                                    include_strand=True, end_sequenced='5prime', print_info=True):
    """ Return all uniquely mappable genomic insertion sites, as dict with chrom or (chrom,strand) keys and pos_array values.

    The mappability is based on the mappability of the flanking region, using flanking_region_length, 
     on either side of the insertion (both sides are counted). 

    The actual positions are the positions of the base BEFORE the insertion, 1-based; 
     the positions are given as numpy arrays (less memory usage than lists), and are sorted. 
    For example if the position 100-120 flanking region is mappable, that means an insertion at 99-100 is mappable, 
     and one at 120-121 is mappable, depending on the orientation: 
      if the flanking region is 5prime, then a +strand 120-121 and a -strand 99-100 insertion is mappable; 
      if 3prime, the strands are inverted. 
    Each cassette orientation gets a separate entry in the position array.
    If include_strand is True, all positions are separated by strand, so the keys in the returned dictionary are
     (chromosome,strand) tuples - so if both a +strand and a -strand insertion at position X of chr1 is mappable, 
      position X will show up once in the ('chr1','+') array and once in ('chr1','-').
     If include_strand is False, the positions are merged, and the return dict keys are just chromosome names,
      so in the same example, position X will just show up twice on the 'chr1' array; 
      if position X is mappable only in one orientation, it'll show up once in the array.
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
    mappable_position_data = collections.defaultdict(list)
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

    # make sure the data is sorted (MAYBE-TODO does it really matter?)
    mappable_position_data = general_utilities.sort_lists_inside_dict(mappable_position_data)
    # Convert the result to a numpy array, it's way faster!  But I have to convert it afterward rather than doing it
    #  as a numpy array to start with, because you can't append to numpy arrays.  MAYBE-TODO or could I do it in-place somehow?  
    mappable_position_data = {key: numpy.array(pos_list) for (key, pos_list) in mappable_position_data.items()}
    return mappable_position_data


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
    # need to start with lists, not arrays, because arrays can't be appended/extended to
    full_mappable_position_data = collections.defaultdict(list)
    # for each set of inputs, run genome_mappable_insertion_sites to generate new data, 
    #  and then merge the new data into full_mappable_position_data
    for flanking_region_length,mappable_slice_pos_dict in zip(flanking_region_lengths, mappable_slice_pos_dicts):
        # get data for 
        curr_mappable_position_data = genome_mappable_insertion_sites(flanking_region_length, mappable_slice_pos_dict, genome_seq, 
                                                                      include_strand, end_sequenced, print_info)
        for chrom,new_pos_list in curr_mappable_position_data.iteritems():
            full_mappable_position_data[chrom].extend(new_pos_list)
    # make sure the data is sorted (MAYBE-TODO does it really matter?) (also converts from defaultdict to dict)
    full_mappable_position_data = general_utilities.sort_lists_inside_dict(full_mappable_position_data)
    # Convert the result to a numpy array, it's way faster!  But I have to convert it afterward rather than doing it
    #  as a numpy array to start with, because you can't append to numpy arrays.  
    full_mappable_position_data = {key: numpy.array(pos_list) for (key, pos_list) in full_mappable_position_data.items()}
    return full_mappable_position_data


############################# Using mappability to finding hotspots/coldspots in real data #############################

### STATISTICAL NOTES (originally from ~/experiments/mutant_pool_screens/1211_positions_Ru-screen1-for-paper/notes.txt)
    # Basic idea - trying to figure out the specific locations of statistically significant hotspots and coldspots:
    # Go over the genome in a sliding window of some size, see how many mutants are there and whether that's statistically significantly different than expected from a random distribution of the same density, INCUDING MAPPABILITY DATA (easiest to do by just comparing to a large number of random datasets, probably).  Figure out if we have any clear hot or cold spots - then we can see if they have anything in common.
    #
    # What statistical test to use for each window?  What we're comparing is the mutant count in the window, X, against one of two things:
    #  A) the mappability % of that window
    #  B) the mutant count in that window for each of the 1000 simulated datasets
    # Which of these things would be the better choice?  Probably A, really, since B is just a random simulation that's entirely based on A.  In either case, REMEMBER TO DO FALSE DISCOVERY RATE CORRECTION after getting the p-values!
    #
    # A) How would I go about getting the p-value for A?  Basically I want to know the probability of getting X mutants in that window randomly, when inserting N mutants in the genome.  That should be mathematically pretty straightforward, right?  All I need is the total mappable size of the genome, and the total mappable size of that window, and I should be able to calculate the actual distribution... That's the binomial distribution - " the discrete probability distribution of the number of successes in a sequence of n independent yes/no experiments, each of which yields success with probability p."  
    # Trying the binomial distribution: the corresponding statistical significance test is the binomial test (https://en.wikipedia.org/wiki/Binomial_test), scipy.stats.binom_test in python.  Say we want the p-value for getting 25 mutants in a 20kb window, out of a total of 12061 mutants, with the probability around 20kb divided by the genome size (really it should be the mappable lengths of that window and the genome, not the total lengths, but let's go with the approximation for now).  The genome size is 113MB, so:
    #     >>> import scipy.stats
    #     >>> scipy.stats.binom_test(25,12061,20000/113000000)
    #     1.3921234115131231e-18
    # (That's pretty significant!  But we'll see how it looks with FDR correction, of course.)
    # How long does this take?  For 50k it takes <1min, and the whole genome in 20bk windows, 113Mb/20kb, is about 5k - not bad.
    #     >>> time.ctime(); _ = [ scipy.stats.binom_test(random.randrange(30),12061,20000/113000000) for x in range(50000)]; time.ctime()
    #     'Sun Dec  2 19:03:25 2012'
    #     'Sun Dec  2 19:04:12 2012'
    # With a reasonably high N (around 1000), the binomial distribution can be modeled by the Poisson distribution (https://en.wikipedia.org/wiki/Poisson_distribution#Derivation_of_Poisson_distribution_.E2.80.94_The_law_of_rare_events).  My N is the number of mutants in the dataset, which is 12k, so I could do that to speed things up - MAYBE-TODO. 
    #
    # B) Alternatively, for option B - I'm trying to test whether the mean of two distributions (real and simulated mutant counts over a window) are the same or different...  So T-test (scipy.stats.ttest_ind, I think), or Welch's t-test (doesn't assume equal variances - apparently present in a newer scipy version), or Mann-Whitney U test (doesn't assume the data is normally distributed).  I'd say it's pretty likely to be randomly distributed, idk about the variances...  Well, I'm only giving a single number for the real dataset, so there's no variance in that - and if anything, the variance of the real dataset may be smaller than of the simulated ones, not larger, so I think we're good with the Student's T-test.
    #
    # How to do FDR correction?  According to the Handbook of Biological Statistics (https://udel.edu/~mcdonald/statmultcomp.html), Benjamini-Hochberg correction is probably what I want.  They describe a procedure, but it's slightly odd, because it doesn't give a p-value (q-value?) for each window, just a yes/no significance result based on the p-value and the desired false discovery rate. 
    # I didn't find any obvious way of doing this directly in python, but there's an R function "p.adjust" (http://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html), which I can use in python with rpy2 (http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python).  Trying that, with just a few test values:
    #  * get the p_values for a few mutant counts per bin, between 0 and 15:
    #     >>> p_values = [scipy.stats.binom_test(x,12061,20000/113000000) for x in (0,0,1,2,5,10,10,25,25)]
    #     >>> p_values
    #     [0.28620628492789102, 0.28620628492789102, 0.73047985928763548, 1.0, 0.065605526425554839, 7.8933016187778668e-05, 7.8933016187778668e-05, 1.3921234115131231e-18, 1.3921234115131231e-18]
    #  * try the FDR adjustment with default options - the p-values increase a bit:
    #     >>> from rpy2.robjects.packages import importr
    #     >>> R_stats = importr('stats')
    #     >>> from rpy2.robjects.vectors import FloatVector
    #     >>> p_adjusted = R_stats.p_adjust(FloatVector(p_values), method = 'BH')
    #     >>> list(p_adjusted)
    #     [0.36797950919300276, 0.36797950919300276, 0.8217898416985899, 1.0, 0.1180899475659987, 0.000177599286422502, 0.000177599286422502, 6.264555351809054e-18, 6.264555351809054e-18]
    #  * same, but using the correct N - I'm reporting 9 random p-values here, but we actually did 50000 tests (50000 windows), not just 9!  The values went down further - good.
    #     >>> p_adjusted2 = R_stats.p_adjust(FloatVector(p_values), method = 'BH', n=50000)
    #     >>> list(p_adjusted2)
    #     [1.0, 1.0, 1.0, 1.0, 1.0, 0.9866627023472333, 0.9866627023472333, 3.4803085287828076e-14, 3.4803085287828076e-14]

def find_hot_cold_spots(dataset, window_size, mappable_positions_20bp, mappable_positions_21bp, window_offset=0, N_samples=None, 
                        chromosome_lengths=None):
    """ Find statistically significant hotspots and coldspots in dataset mutant positions, based on genome mappability. 

    Divides the genome into window_size-sized windows, and for each of them compares 
     the number of mutants in that window in dataset to the number expected using the BINOMIAL DISTRIBUTION, 
      given the total number of mutants in dataset and the genome mappability information
      (the probability of a mutant landing in the window is the number of mappable positions in the window 
       divided by the number of mappable positions in the entire genome - these numbers will be different for 20bp 
       and 21bp mappability, so calculate for both, and get a weighted average based on how many mutants in dataset
       have only 20bp flanking region sequences.)
      FALSE DISCOVERY RATE CORRECTION using the Benjamini-Hochberg method is used on all p-values, 
       with the total number of windows tested used as N_samples, unless another N_samples is provided
       (which you may or may not want to do if doing multiple find_hot_cold_spots runs 
        with different window sizes/offsets over the same dataset: you'll be under-correcting if you don't, 
        but since the results for overlapping windows are positively correlated, you'll be over-correcting if you do!)

    Window_offset defines at which position of each chromosome to start: if a chromosome is 500 long and window_size is 200, 
     - with window_offset 0, the windows checked will be 1-200 and 201-400;
     - with window_offset 100, the windows checked will be 101-300 and 301-500.
    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default chlamy file will be used) - need the lengths because the dataset positions and mappability data
     don't give the end position of each chromosome.

    Return a 4-tuple of FDR-adjusted p-values, raw p-values, sides (1 if we had more mutants than expected, -1 if fewer),
      and raw mutant counts, with each of these items being a chromosome:list_of_values_per_window dictionary 
     the caller/recipient must take care of keeping track of the window size and offset in order for this data to be meaningful.

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance. 

    The two mappable_positions_* arguments should be (chrom,strand):position_list dictionaries 
     giving all the mappable positions, as generated by genome_mappable_insertion_sites with either include_strand value. 
    """
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
    fraction_20bp = get_20bp_fraction(dataset)
    total_N_mutants = len(dataset)
    R_stats = importr('stats')
    mappable_position_data = {20: mappable_positions_20bp, 21: mappable_positions_21bp}
    # get total mappable genome lengths (for 20 and 21bp cases separately)
    genome_mappable_lengths = {flank_len: sum(len(pos_list) for pos_list in mappable_positions.values())
                               for (flank_len,mappable_positions) in mappable_position_data.items()}
    # get mappable lengths per window, for 20 and 21bp cases separately, with the right window size and offset; 
    #  don't include special_last_bin to avoid complications (so the last part of each chromosome that doesn't 
    #   fit the window size will be ignored - same as the first part before offset)
    #   (chromosomes smaller than window_size+offset will be ignored entirely!)
    window_kwargs = {'bin_size':window_size, 'chromosome_lengths':chromosome_lengths, 
                     'first_bin_offset': window_offset, 'special_last_bin':False}
    window_mappable_length_lists = {flank_len: get_histogram_data_from_positions(mappable_positions, **window_kwargs) 
                                    for (flank_len,mappable_positions) in mappable_position_data.items()}
    # similarly, get mutant counts per window
    window_mutant_count_lists = get_histogram_data_from_positions(get_mutant_positions_from_dataset(dataset), **window_kwargs)
    total_N_windows = sum(len(len_list) for len_list in window_mutant_count_lists.values())
    if N_samples is None:   N_samples = total_N_windows
    window_raw_pvalues = {}
    window_adjusted_pvalues = {}
    window_which_side_values = {}
    # average mutants per window should only count the mutants covered by any of the windows 
    #  (exclude ones skipped due to non-zero offset at chromosome start or a non-full window at chromosome end)
    average_mutants_per_window = sum(sum(counts) for counts in window_mutant_count_lists.values()) / total_N_windows
    print "%s window (offset %s) - average %.2g mutants per window "%(basic_seq_utilities.format_base_distance(window_size, False), 
                                          basic_seq_utilities.format_base_distance(window_offset, False), average_mutants_per_window)
    for chrom, window_mutant_count_list in window_mutant_count_lists.items():
        # get the probability of a mutant landing in the window, 
        #  given the mappable lengths of the window and the total genome for 20 and 21bp cases, 
        #  and the fraction of mutants that is 20bp-only.
        # note that all the values in the histogram dicts are numpy arrays, not lists, so they can be operated on directly
        window_probabilities_20 = window_mappable_length_lists[20][chrom] / genome_mappable_lengths[20]
        window_probabilities_21 = window_mappable_length_lists[21][chrom] / genome_mappable_lengths[21]
        window_probabilities = window_probabilities_20 * fraction_20bp + window_probabilities_21 * (1-fraction_20bp)
        window_which_side_values[chrom] = numpy.array([-1 if m < average_mutants_per_window else 1 
                                                       for m in window_mutant_count_list])
        # calculate p-value for each window, using the binomial distribution 
        #  based on the window probability and window mutant number; convert to numpy array again
        raw_pvalues = [scipy.stats.binom_test(x=N_mutants_in_window, n=total_N_mutants, p=window_probability) 
                       for (N_mutants_in_window, window_probability) in zip(window_mutant_count_list, window_probabilities)]
        window_raw_pvalues[chrom] = numpy.array(raw_pvalues)
        # adjust the p-values for FDR, using the total N_samples (rather than just the number of samples in this chromosome)
        #  (if there are no windows on a given chromosome (it was too short for the size+offset), 
        #   just set the lists to empty and skip rather than running statistics on empty lists, which can give warnings)
        if len(window_mutant_count_list):
            window_adjusted_pvalues[chrom] = numpy.array(R_stats.p_adjust(FloatVector(window_raw_pvalues[chrom]), 
                                                                          method='BH', n=N_samples))
        else:
            window_adjusted_pvalues[chrom] = numpy.array([])
    return window_adjusted_pvalues, window_raw_pvalues, window_which_side_values, window_mutant_count_lists
    # TODO unit-test this? How?


def get_hot_cold_spot_list(pvalue_data_window_size_offset_dict, side_data_window_size_offset_dict, 
                           pval_cutoff=0.05, print_info=True):
    """ Transform p-value and side info into a list of (chrom, start, end, pvalue, kind, window_offset) tuples with pvalue<=cutoff.

    Both arguments should be window_size:(window_offset:(chromosome:list_of_window_values))) triple nested dictionaries.
      - the p-values should be 0-1, and you should probably get FDR correction done on them first; 
      - the sides should be -1 for coldspots and 1 for hotspots. 
    Pval_cutoff should be a number between 0 and 1. 

    In the output, kind is 'hotspot' or 'coldspot'. 
    Sort the output data by chromosome/position;  If print_info, print the output data.
    """
    hc_data_list = []
    for window_size, pvalue_data_window_offset_dict in pvalue_data_window_size_offset_dict.items():
        for window_offset, pvalue_data in pvalue_data_window_offset_dict.items():
            side_data = side_data_window_size_offset_dict[window_size][window_offset]
            for chrom, pvalues in pvalue_data.items():
                sides = side_data[chrom]
                color_values = []
                for N,(pvalue,side) in enumerate(zip(pvalues,sides)):
                    if pvalue <= pval_cutoff:
                        kind = 'hotspot' if side>0 else 'coldspot'
                        start_pos = N*window_size + window_offset
                        end_pos = (N+1)*window_size + window_offset
                        hc_data_list.append((chrom, start_pos, end_pos, pvalue, kind, window_offset))
    hc_data_list.sort(key=lambda (c,s,e,p,k,o): (basic_seq_utilities.chromosome_sort_key(c), s,e,p,k,o))
    if print_info:
        format_bp = lambda x: basic_seq_utilities.format_base_distance(x, False)    # the False is to not approximate
        for (chrom, start_pos, end_pos, pvalue, kind, offset) in hc_data_list:
            print "%s %s-%s (window size %s) - %s, %.3g adjusted p-value"%(chrom, format_bp(start_pos), format_bp(end_pos), 
                                                                           format_bp(end_pos-start_pos), kind, pvalue)
    return hc_data_list
    # TODO should unit-test this!


################################# Simulating random size-N dataset, various options #######################################

def weighted_random_choice(value_weight_list):
    """ Given a list of (value, weight) tuples, pick a random value, with probability=weight for each value. """
    # also see http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/ - using weighted_choice_sub
    #  MAYBE-TODO make this faster by using the WeightedRandomGenerator class from the same site?
    #  also see http://docs.python.org/3/library/random.html - similar, but for python3 only, (2.7 is missing itertools.accumulate)
    values, weights = zip(*value_weight_list)
    rnd = random.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return values[i]
    # LATER-TODO move this to general_utilities?


def simulate_dataset_from_mappability(N_mutants, fraction_20bp, mappable_positions_20bp, mappable_positions_21bp, 
                                      all_chromosomes=None, include_strand=False, fraction_plus_strand=0.5):
    """ Return position info from simulating inserting N mutants into the genome randomly, taking mappability into account. 
    
    Return a (chrom,strand):position_list dictionary if include_strand, or else just a chrom:position_list one...
     (the latter is useful for plotting or gap size analysis, since both of those ignore strand.

    Fraction_20bp should be the fraction of mutants that have 20bp vs 21bp flanking regions 
     (since the mappability data is slightly different for those two cases).
    Fraction_plus_strand is the fraction of mutants on the +strand (vs the -strand); 
     note that this is STILL RELEVANT even if include_strand is False, because mappability differs slightly between strands!

    All_chromosomes should give a list of chromosomes to include (e.g. if trying to simulate a dataset without chloroplast/mito
     genomes, don't put those on a list!); if None, all chromosomes in mappable_positions_20bp will be used.

    The two mappable_positions_* arguments should be (chrom,strand):position_list dictionaries given all the mappable positions,
     as generated by genome_mappable_insertion_sites with include_strand==True.
    """
    # MAYBE-TODO add a way of copying N_mutants and fraction_* from an existing dataset instead of putting them in separately?  Could do that by putting it here as an argument, or with a help function or something
    if all_chromosomes is None:
        all_chromosomes = set(chrom for chrom,s in mappable_positions_20bp) | set(chrom for chrom,s in mappable_positions_21bp)
    mappable_position_data = {20: mappable_positions_20bp, 21: mappable_positions_21bp}
    # chromosome mappable lengths - a list of (chrom,mappable_lengths) tuples for each (flank_length,chrom,strand) combination, 
    #  to use as a list of values and weights when randomly choosing a chromosome - only include ones in all_chromosomes!
    chrom_mappable_len = collections.defaultdict(list)
    for flank_len,mappable_positions in mappable_position_data.iteritems():
        for (chrom,strand),pos_list in mappable_positions.iteritems():
            if chrom in all_chromosomes:
                chrom_mappable_len[flank_len,strand].append((chrom, len(pos_list)))
    simulated_positions = collections.defaultdict(list)
    for _ in range(N_mutants):
        # first choose a length (20 or 21bp) based on Fraction_20bp, and a strand based on fraction_plus_strand.
        #   random.random() gives a value x in the interval [0, 1) (including 0, excluding 1).
        #   if fraction_* is 0, we want to never pick that option, so the right test is random<fraction.
        if random.random() < fraction_20bp: flank_length = 20
        else:                               flank_length = 21
        if random.random() < fraction_plus_strand: strand = '+'
        else:                                      strand = '-'
        # next choose a chromosome, with a probability corresponding to the mappable length of each chromosome
        #  (given the flank_length and the strand)
        chrom = weighted_random_choice(chrom_mappable_len[(flank_length,strand)])
        # next choose a position from the list of mappable positions on the chromosome
        pos = random.choice(mappable_position_data[flank_length][(chrom,strand)])
        # save the data
        if include_strand:  simulated_positions[(chrom,strand)].append(pos)
        else:               simulated_positions[chrom].append(pos)
    # convert output to numpy arrays for speed and memory efficiency
    for key,val_list in simulated_positions.iteritems():
        simulated_positions[key] = numpy.array(val_list)
    return simulated_positions
    # TODO how would I test this sensibly??  Complicated... But I did run it, visualize the results, and compare to a real dataset, and it looked all right.  (../1211_positions_Ru-screen1-for-paper/positions_over_chromosomes/*simulated*)


### simulate dataset with N randomly positioned mutants, matching the gap-size distribution of the real dataset, and taking into account mappable/unmappable positions
# LATER-TODO
# MAYBE-TODO would it be possible to also match the real hotspots and coldspots of the real dataset?


def find_genes_for_mutants(self, genefile, detailed_features=True, N_run_groups=3, verbosity_level=1):
    """ To each mutant in the dataset, add the gene it's in (look up gene positions for each mutant using genefile).

    If detailed_features is True, also look up whether the mutant is in an exon/intron/UTR (NOT IMPLEMENTED); 
    Read the file in N_run_groups passes to avoid using up too much memory/CPU.

    Based on mutant_analysis_classes.Insertional_mutant_pool_dataset.find_genes_for_mutants,
     but acts on simple position tuples instead of mutant objects, and just returns the list of genes.

    

    """ 
    gene_list = []
    ### go over all mutants on each chromosome, figure out which gene they're in (if any), keep track of totals
    # keep track of all the mutant and reference chromosomes to catch chromosomes that are absent in reference
    with open(genefile) as GENEFILE:
        for chromosome_record in GFF.parse(GENEFILE, limit_info=genefile_parsing_limits):
            self.total_genes_in_genome += len(chromosome_record.features)
            # TODO TODO TODO implement this bit!
            for mutant in mutants_by_chromosome[chromosome_record.id]:
                gene_ID, orientation, feature = find_gene_by_pos(mutant.position, chromosome_record, 
                                                                 detailed_features, quiet=(verbosity_level==0))
                mutant.gene, mutant.orientation, mutant.gene_feature = gene_ID, orientation, feature

    # for mutants in chromosomes that weren't listed in the genefile, use special values
    for chromosome in set(mutants_by_chromosome.keys())-set(all_reference_chromosomes):
        if not is_cassette_chromosome(chromosome):
            print 'Warning: chromosome "%s" not found in genefile data!'%(chromosome)
        for mutant in mutants_by_chromosome[chromosome]:
            mutant.gene,mutant.orientation,mutant.gene_feature = SPECIAL_GENE_CODES.chromosome_not_in_reference,'-','-'


def dataset_object_from_simulated(simulated_dataset, get_gene_info=True, 
                  gene_info_file=os.path.expanduser("~/experiments/reference_data/chlamy_annotation/Creinhardtii_169_gene.gff3")):
    """ Convert a simulated dataset position dictionary into an Insertional_mutant_pool_dataset object, with optional gene info.
    
    Output will be a mutant_analysis_classes.Insertional_mutant_pool_dataset object, with all readcounts set to 1.

    If get_gene_info is True, get the gene information for the positions (from gene_info_file if given, otherwise the default).

    Input should be a (chrom,strand):position_list dictionary, like from simulate_dataset_from_mappability with include_strand True.

    Warning: can take a LOT of time/memory if the simulated dataset is large!  May be best to avoid if possible.
    """
    new_dataset = mutant_analysis_classes.Insertional_mutant_pool_dataset()
    # Go over all the simulated read positions, add them to dataset
    for (chrom,strand),pos_list in simulated_dataset.items():
        for pos in pos_list:
            read_count = 1
            # make a simple position - just always specify position_before, regardless of strand
            position = mutant_analysis_classes.Insertion_position(chrom, strand, position_before=pos, immutable=True)
            # grab the right mutant based on the position, and add the reads to it; 
            curr_mutant = new_dataset.get_mutant(position)
            # just add 1 readcount, and don't bother with the sequences etc
            curr_mutant.add_counts(1,0,1)
    # Add gene info to dataset if desired
    if get_gene_info:
        new_dataset.find_genes_for_mutants(gene_info_file, detailed_features=True)
    return new_dataset


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
    # TODO TODO TODO modify this to just extract a list of GENES and work on that!
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
            formatted_data = collections.defaultdict(list)
            for chr_data in raw_data.split(', '):
                chrom, pos_data = chr_data.strip().split(': ')
                # pos_data can be ints (1 25 301) or ints with strand info (1- 25+ 301+)
                for x in pos_data.split(' '):
                    if x[-1] in '+-':   formatted_data[chrom,x[-1]].append(int(x[:-1]))
                    else:               formatted_data[chrom].append(int(x))
            for val in formatted_data.values():
                val.sort()
            return dict(formatted_data)
        def _compare_dicts(numpy_array_dict, list_dict):
            """ Compare a key:val_numpy_array dict to a key:val_list dict, make sure they match. """
            assert numpy_array_dict.keys() == list_dict.keys()
            for key,val in numpy_array_dict.iteritems():
                assert list(val) == list_dict[key]
        def _test_all(slice_len, genome, raw_slice_data, raw_pos_data_5prime):
            """ Test genome_mappable_slices and all variants of genome_mappable_insertion_sites* against expected output. """
            # get the expected output data from the simplified string formats
            slice_data = _read_raw_data(raw_slice_data)
            pos_data_5prime = _read_raw_data(raw_pos_data_5prime)
            # from pos_data_5prime, make pos_data_3prime (just switch all strands) and pos_data_no_strand (remove strand info)
            pos_data_3prime, pos_data_no_strand = {}, collections.defaultdict(list)
            for (chrom,strand),pos_list in pos_data_5prime.items():
                pos_data_3prime[chrom, '+' if strand=='-' else '-'] = pos_list
                pos_data_no_strand[chrom] += pos_list
                pos_data_no_strand[chrom].sort()
            # check genome_mappable_slices output (and save it for later)
            new_slice_data = genome_mappable_slices(slice_len, genome, False)
            _compare_dicts(new_slice_data, slice_data)
            # now try running genome_mappable_insertion_sites with both the raw slice_len/genome data, and the new_slice_data; 
            for include_strand,end,pos_data in ((True,'5prime',pos_data_5prime), (True,'3prime',pos_data_3prime), 
                                                (False,'5prime',pos_data_no_strand), (False,'3prime',pos_data_no_strand)):
                    args = (include_strand,end, False)
                    _compare_dicts(genome_mappable_insertion_sites(slice_len, slice_data, None, *args), pos_data)
                    _compare_dicts(genome_mappable_insertion_sites(slice_len, None, genome, *args), pos_data)
                    _compare_dicts(genome_mappable_insertion_sites_multi([slice_len], None, genome, *args), pos_data)
                    for extra in ([], [{}], [{}, {}, {}]):
                        _compare_dicts(genome_mappable_insertion_sites_multi([slice_len]*(len(extra)+1), [slice_data]+extra, 
                                                                                    None, *args), pos_data)
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
        _compare_dicts(genome_mappable_insertion_sites_multi([1,2], None, curr_genome, False, '5prime', False), {'a':[0,1,1,2,2,3]})
        slices_1 = genome_mappable_slices(1, curr_genome, False)
        slices_2 = genome_mappable_slices(2, curr_genome, False)
        for fl_both,slices_both in (([1,2],[slices_1,slices_2]), ([2,1],[slices_2,slices_1])):
            _compare_dicts(genome_mappable_insertion_sites_multi(fl_both, slices_both, None, False, '5prime', False), 
                           {'a':[0,1,1,2,2,3]})

        
    # LATER-TODO add unit-tests for other stuff!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
