#! /usr/bin/env python

"""
Various utilities for analysis of mutant flanking regions (in progress, may not be very clean).
 -- Weronika Patena, 2013
"""

# standard library
from __future__ import division
import os, sys
import unittest
import itertools
from collections import defaultdict
import random
# other packages
import scipy.stats
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties
# my modules
import mutant_analysis_classes
from basic_seq_utilities import get_all_seq_length, base_count_dict, base_fraction_dict, base_fraction_dict_from_count_dict, base_fractions_from_GC_content, NORMAL_DNA_BASES, reverse_complement, check_seq_against_pattern, write_fasta_line
import general_utilities
import statistics_utilities


def flanking_region_from_pos(position_before_insertion, chromosome, strand, genome, flanksize=200, padding_char='.'):
    """ Given a position (pos,chromosome,strand), and the genome sequence, return position flanking region

    Given a position, chromosome and strand for an insertion, use genome (a chrom_name:seq dict) to figure out the 
     flanksize-length flanking sequences on both sides (padded with padding_char if the end of the chromosome is too close), 
     reverse-complement if needed (if strand=='-') to get it in the same orientation as the insertion, and return as string.
    """
    chromosome_seq = genome[chromosome]
    # grab positions to get the sequence flanksize bp on both sides of the insertion;
    #  pad it with . if there's a chromosome start/end closer than that.
    flank_start = position_before_insertion - flanksize
    flank_end = position_before_insertion + flanksize
    start_padding = padding_char * max(0, -flank_start)
    flank_start = max(flank_start, 0)
    end_padding = padding_char * max(0, flank_end - len(chromosome_seq))
    flank_end = min(flank_end, len(chromosome_seq))
    # grab the actual flanking sequence, add padding; 
    #  reverse-complement it if mutant is -strand, to keep them in the same orientation as the cassette insertion
    full_flanking_seq = chromosome_seq[flank_start:flank_end]
    full_flanking_seq = start_padding + full_flanking_seq + end_padding
    if strand == '-':
        full_flanking_seq = reverse_complement(full_flanking_seq)
    return full_flanking_seq
    # TODO unit-test this!


def grab_flanking_regions_from_mutantfile(mutant_dataset_infile, genome, flanksize=200, padding_char='.', 
                                      min_readcount=0, chromosome_check_function=lambda x: True, ignore_both_strand_mutants=False):
    """ Return (flanking_seq,readcount) with both-side genomic flanking sequences for insertional mutants in mutant_dataset_infile.

    Grab all the insertion positions from mutant_dataset_infile (pickled mutant_analysis_classes.Insertional_mutant_dataset object), 
     use genome (a chrom_name:seq dict) to figure out the flanksize-length flanking sequences on both sides
      (padded with padding_char if the end of the chromosome is too close), reverse-complement if needed (if strand=='-') 
      to get it in the same orientation as the insertion.

    Filter the mutants: 
     - by readcount - ignore mutants with total readcount below min_readcount=0
     - by chromosome - ignore mutants in chromosomes for which chromosome_check_function returns False
     - by strand - both-strand (merged tandem) mutants will be ignored if ignore_both_strand_mutants is True, 
                otherwise ValueError will be raised; ValueError will be raised for other unexpected strand values.

    For all remaining mutants, append (flanking region seq, total_readcount) to output list.
    """
    dataset = mutant_analysis_classes.read_mutant_file(mutant_dataset_infile)
    flanking_region_count_list = []
    for mutant in sorted(dataset, key = lambda m: m.position):
        # filter out mutants with wrong readcounts or in wrong chromosomes
        if not chromosome_check_function(mutant.position.chromosome):           continue
        if mutant.total_read_count < min_readcount:                             continue
        # filter out both-stranded mutants if desired; 
        if mutant.position.strand not in '+-':  
            if mutant.position.strand=='both' and ignore_both_strand_mutants:   continue
            else:                                   raise ValueError("Unexpected mutant strand! %s"%mutant.position)
        # grab mutant position/chromosome
        position_before_insertion = mutant.position.min_position
        # ignore cassette tandems (i.e. insertions that map to start or end of cassette)
        if mutant_analysis_classes.is_cassette_chromosome(mutant.position.chromosome):
            if position_before_insertion in [0, len(genome[mutant.position.chromosome])]:  continue
        # grab the actual flanking sequence, with padding, correct orientation etc
        full_flanking_seq = flanking_region_from_pos(position_before_insertion, mutant.position.chromosome, 
                                                     mutant.position.strand, genome, flanksize, padding_char)
        # append the sequence and readcount to output data
        flanking_region_count_list.append((full_flanking_seq, mutant.total_read_count))
    return flanking_region_count_list


def grab_flanking_regions_from_pos_dict(insertion_position_dict, genome, flanksize=200, padding_char='.', 
                                        chromosome_check_function=lambda x: True, ignore_both_strand_mutants=False):
    """ Same as grab_flanking_regions_from_mutantfile, but takes input as a (chrom,strand):pos_list dictionary instead, 
    and assumes all readcounts to be 1.
    """
    flanking_region_count_list = []
    for (chromosome,strand),pos_list in insertion_position_dict.items():
        for position_before_insertion in pos_list:
            # filter out positions in wrong chromosomes; filter out both-stranded positions if desired
            if not chromosome_check_function(chromosome):           continue
            if strand not in '+-':  
                if strand=='both' and ignore_both_strand_mutants:   continue
                else:                                               raise ValueError("Unexpected strand! %s"%strand)
            # ignore cassette tandems (i.e. insertions that map to start or end of cassette)
            if mutant_analysis_classes.is_cassette_chromosome(chromosome):
                if position_before_insertion in [0, len(genome[chromosome])]:  continue
            # grab the actual flanking sequence, with padding, correct orientation etc
            full_flanking_seq = flanking_region_from_pos(position_before_insertion, chromosome, strand, 
                                                         genome, flanksize, padding_char)
            # append the sequence and readcount to output data
            flanking_region_count_list.append((full_flanking_seq, 1))
    return flanking_region_count_list
    # TODO unit-test this too?  Could probably just convert the mutant data from the grab_flanking_regions_from_mutantfile unit-test and use that.


def grab_flanking_region_base_counts_from_pos_dict(insertion_position_dict, genome, flanksize=200, 
                                                   chromosome_check_function=lambda x: True, ignore_both_strand_mutants=False):
    """ Same as base_count_dict(grab_flanking_regions_from_mutantfile(*args)), but saves memory by not keeping all the sequences.

    Basically instead of making a full dataset of flanking regions with grab_flanking_regions_from_mutantfile (which can be BIG)
     and then converting those to a base-count dict with base_count_dict, just go over each position in insertion_position_dict,
     grab that flanking region, add it to the current base-count dict, and go on to the next one, without saving.
    Assumes all readcounts to be 1.
    """
    # initialize the base-count lists to the right length, and fill it out by going over all the seqs
    base_count_dict = {base: [0 for _ in range(flanksize*2)] for base in NORMAL_DNA_BASES}
    # for each position, grab the flanking region and add it to the base_count_dict
    for (chromosome,strand),pos_list in insertion_position_dict.items():
        for position_before_insertion in pos_list:
            # filter out positions in wrong chromosomes; filter out both-stranded positions if desired
            if not chromosome_check_function(chromosome):           continue
            if strand not in '+-':  
                if strand=='both' and ignore_both_strand_mutants:   continue
                else:                                               raise ValueError("Unexpected strand! %s"%strand)
            # ignore cassette tandems (i.e. insertions that map to start or end of cassette)
            if mutant_analysis_classes.is_cassette_chromosome(chromosome):
                if position_before_insertion in [0, len(genome[chromosome])]:  continue
            # grab the actual flanking sequence, with padding, correct orientation etc
            full_flanking_seq = flanking_region_from_pos(position_before_insertion, chromosome, strand, genome, flanksize)
            # add base-counts from full_flanking_seq to base_count_dict
            for position, base in enumerate(full_flanking_seq.upper()):
                try:
                    base_count_dict[base][position] += 1
                except KeyError:
                    pass
                    # MAYBE-TODO add an option to NOT ignore bases that aren't in NORMAL_DNA_BASES?
    return base_count_dict
    # TODO unit-test this too?


def grab_flanking_region_motif_counts_from_pos_dict(insertion_position_dict, genome, flanksize=2, 
                                                    chromosome_check_function=lambda x: True, ignore_both_strand_mutants=False):
    """ Get a flanking_seq:count dictionary for the flanking seqs for insertion_position_dict ((chrom,strand):pos_list dictionary).

    Only really makes sense for small flanksizes - otherwise the total number of possible motifs will be huge (4^(flanksize*2).
    Assumes all readcounts to be 1.
    """
    # initialize the motif-count lists to the right length, and fill it out by going over all the seqs
    motif_count_dict = {''.join(four_bases): 0 for four_bases 
                        in itertools.product(NORMAL_DNA_BASES,NORMAL_DNA_BASES,NORMAL_DNA_BASES,NORMAL_DNA_BASES)}
    # for each position, grab the flanking region and add it to the motif_count_dict
    for (chromosome,strand),pos_list in insertion_position_dict.items():
        for position_before_insertion in pos_list:
            # filter out positions in wrong chromosomes; filter out both-stranded positions if desired
            if not chromosome_check_function(chromosome):           continue
            if strand not in '+-':  
                if strand=='both' and ignore_both_strand_mutants:   continue
                else:                                               raise ValueError("Unexpected strand! %s"%strand)
            # ignore cassette tandems (i.e. insertions that map to start or end of cassette)
            if mutant_analysis_classes.is_cassette_chromosome(chromosome):
                if position_before_insertion in [0, len(genome[chromosome])]:  continue
            # grab the actual flanking sequence, with padding, correct orientation etc
            full_flanking_seq = flanking_region_from_pos(position_before_insertion, chromosome, strand, genome, flanksize)
            # add motif-count of full_flanking_seq to motif_count_dict
            try:                motif_count_dict[full_flanking_seq] += 1
            except KeyError:    pass
            # MAYBE-TODO add an option to NOT ignore motifs that aren't in NORMAL_DNA_BASES?
    return motif_count_dict
    # TODO unit-test this too?
    # TODO should probably just merge this with grab_flanking_region_base_counts_from_pos_dict for the future?  Except I'd usually want to run them with different flanksizes, so maybe not.  May want to refactor or something - give different flanksizes for the two functionalities but in a single function?  Since most of the work is probably grabbing the flanking region...


def filter_flanking_regions_by_pattern(flanking_region_count_list, pattern, either_orientation=True, 
                                       print_info=True, category=None, meaning_of_seqs='positions', meaning_of_counts='counts'):
    """ Return separate lists of flanking regions that do and don't match given sequence pattern.
    
    flanking_region_count_list should be a list of (flanking_region, count) pairs (like from grab_flanking_regions_from_mutantfile); 
     the two return values (flanking regions that match and don't match the pattern) are the same format.

    The pattern should be a sequence string (allowed letters are ACTGN). It'ss considered to be centered around the cut site; 
     the flanking regions likewise.  E.g. if pattern is GNAN, a fl.region of GCAC or TTGCACTT would match, but TTTTGCAC would not. 
    If either_orientation is True, each flanking region will be tried against the pattern in both the forward and the reverse
     orientation, and the returned flanking region will be in the orientation that matched - e.g. if pattern is GNAN, 
     a flanking region of either TTGCACTT or TTCTCCTT would match (forward and rev-compl respectively), 
      and the latter would be returned as rev-compl, AAGGAGAA.

    If print_info is True, some information will be printed about what number/percentage matched and didn't: it'll be given two ways:
     - by flanking region, counting each once, if meaning_of_seqs is not None, and meaning_of_seqs will be used as the description
     - by count, if some counts are not 1 and meaning_of_counts is not None, and meaning_of_counts will be used as the description.
    """
    if not flanking_region_count_list:  return []
    flanking_region_length = get_all_seq_length(zip(*flanking_region_count_list)[0])
    if flanking_region_length % 2:              raise ValueError("Flanking region length must be an even number!")
    if len(pattern) % 2:                        raise ValueError("Pattern length must be an even number!")
    if len(pattern) > flanking_region_length:   raise ValueError("Pattern cannot be longer than flanking regions!")
    # pad the pattern to match the flanking region length
    orig_pattern = pattern
    if len(pattern) < flanking_region_length:
        padding_len = int((flanking_region_length - len(pattern)) / 2)
        pattern = 'N'*padding_len + pattern + 'N'*padding_len
    # go over all the flanking regions: 
    flanking_region_count_list_match, flanking_region_count_list_nomatch = [], []
    for (flanking_region,count) in flanking_region_count_list:
        # if the flanking region is padded with .'s, change them to N's to make check_seq_against_pattern take it
        flanking_region = flanking_region.replace('.', 'N')
        # if we're looking at both orientations, then first randomize the orientation to avoid bias
        if either_orientation and random.random() < 0.5:
            flanking_region = reverse_complement(flanking_region)
        # if it matches the pattern, save it as a match and go to the next one 
        if check_seq_against_pattern(flanking_region, pattern):
            flanking_region_count_list_match.append((flanking_region, count))
            continue
        # or if its rev-compl matches the pattern and either_orientation is True, save it as a match and go on to the next one; 
        if either_orientation:
            flanking_region = reverse_complement(flanking_region)
            if check_seq_against_pattern(flanking_region, pattern):
                flanking_region_count_list_match.append((flanking_region, count))
                continue
        # if it didn't match anywhere, save it as a no-match.
        flanking_region_count_list_nomatch.append((flanking_region, count))
    if print_info:
        if meaning_of_seqs is None and meaning_of_counts is None:
            raise ValueError("To get info printed, at least one of meaning_of_seqs/meaning_of_counts must be not None!")
        print_data = "%smatched %s:  "%('' if category is None else category+' ', orig_pattern)
        if meaning_of_seqs is not None:
            positions_matched, positions_unmatched = len(flanking_region_count_list_match), len(flanking_region_count_list_nomatch)
            positions_all = positions_matched+positions_unmatched
            print_data += "%s, unmatched %s/%s"%( 
                            general_utilities.value_and_percentages(positions_matched,[positions_all],insert_word=meaning_of_seqs),
                            positions_unmatched, positions_all)
        if meaning_of_counts is not None:
            counts_matched, counts_unmatched = [sum(zip(*data)[1]) 
                                                for data in (flanking_region_count_list_match,flanking_region_count_list_nomatch)]
            counts_all = counts_matched+counts_unmatched
            print_data += ";  %s, unmatched %s/%s."%(
                                general_utilities.value_and_percentages(counts_matched,[counts_all],insert_word=meaning_of_counts),
                                counts_unmatched, counts_all)
        print print_data
    return flanking_region_count_list_match, flanking_region_count_list_nomatch


def print_flanking_regions_to_fasta(flanking_region_count_list, outfile, convert_counts=lambda x: x):
    """ Given a (seq,count) list, make fasta file with the each seq present convert_counts(count) times. """
    with open(outfile, 'w') as OUTFILE:
        for N,(seq,count) in enumerate(flanking_region_count_list):
            seqname = '%s (%s reads)' % (N, count)
            for _ in range(convert_counts(count)):
                write_fasta_line(seqname, seq, OUTFILE)


def _relative_position_vs_cut(pos, reference_pos):
    """ Convert zero-based position into relative position around a cut site before reference_pos, with no 0.

    If refpos is 2, we assume the true reference cut position is between 1 and 2, so positions 1 and 2 become -1 and 1, 
     0 and 3 become -2 and 2, etc.
    """
    if pos < reference_pos:  return int(pos - reference_pos)
    else:                    return int(pos - reference_pos + 1)


def print_base_count_fraction_for_dist(flanking_region_count_list, distance_from_middle, convert_counts=lambda x: x, 
                                      ignore_bases_pattern=None, average_both_sides=False):
    """ Print the base counts/fractions for the N bases around the middle, for all the flanking regions, weighed by converted count.

    Flanking_region_count_list should be a (seq,count) list; each sequence will be weighed as convert_counts(count). 

    Ignore_bases_pattern should be None to count all the bases, or a string giving a base pattern (centered around the middle), 
     in which case only the bases with an N will be counted - so for instnace if ignore_bases_pattern is CANN, bases -2 and -1 
      will be ignored (presumed to have been filtered through a pattern that requires them to be CA, although this isn't checked), 
      and only data for bases 1 and 2 (as well as bases before -2 and after 2, depending on distance_from_middle) will be given.

    If average_both_sides is True, bases from the two sides around the middle will be averaged together (after reverse-complementing
     one side, of course): CA|GG will be converted to |GG and |TG (rev-compl of CA|) and the GG and TG treated together for 
     calculating base frequencies/counts.
    """
    flanking_region_length = get_all_seq_length(zip(*flanking_region_count_list)[0])
    if flanking_region_length % 2:              raise ValueError("Flanking region length must be an even number!")
    flank_length = int(flanking_region_length/2)
    # grab only the flanking region length we're actually interested in, and convert the counts
    local_flanking_region_length = 2*distance_from_middle
    local_flanking_region_count_list = [(seq[flank_length-distance_from_middle:flank_length+distance_from_middle], 
                                         convert_counts(count)) for (seq,count) in flanking_region_count_list]
    # apply ignore_bases_pattern by changing the relevant bases to N
    #  note that I'm doing it this way, instead of just skipping these positions in the final output, 
    #  because that wouldn't work if average_both_sides is True: if ignore_bases_pattern is ANAN, 
    #   the ignore pattern isn't symmetrical around the middle, so half the position 1 bases will be ignored (first N in ANAN), 
    #    and half the position 2 bases will be ignored (second N in ANAN) - there will be no single position with all bases ignored.
    if ignore_bases_pattern is not None:
        if len(ignore_bases_pattern) % 2:       raise ValueError("Ignore_bases_pattern length must be an even number!")
        length_diff = int((len(ignore_bases_pattern)-local_flanking_region_length)/2)
        if length_diff>0:   raise ValueError("Ignore_bases_pattern is longer than 2*distance_from_middle - probably error!")
        if length_diff<0:   ignore_bases_pattern = 'N'*length_diff + ignore_bases_pattern + 'N'*length_diff
        def mask_seq(seq, mask_pattern):
            return ''.join([(base if if_mask=='N' else 'N') for (base,if_mask) in zip(seq, mask_pattern)])
        local_flanking_region_count_list = [(mask_seq(seq, ignore_bases_pattern), count) 
                                            for (seq,count) in local_flanking_region_count_list]
    # if average_both_sides, make a new local_flanking_region_count_list that has each half of each sequence separately
    if average_both_sides:
        new_flanking_region_count_list = []
        for flanking_region,count in local_flanking_region_count_list:
            first_half = reverse_complement(flanking_region[:distance_from_middle])
            second_half = flanking_region[distance_from_middle:]
            new_flanking_region_count_list.extend([(first_half,count), (second_half,count)])
        local_flanking_region_count_list = new_flanking_region_count_list
        local_flanking_region_length = int(local_flanking_region_length/2)
    base_count_list_dict = base_count_dict(local_flanking_region_count_list)
    base_fraction_list_dict = base_fraction_dict(local_flanking_region_count_list)
    # for each position in the final flanking regions, give the base fraction/count; 
    #  ignore positions in which there were no non-N bases.
    all_lines = ''
    for position in range(local_flanking_region_length):
        if sum(base_count_list_dict[base][position] for base in NORMAL_DNA_BASES):
            data = ['%s %.0f%% (%s)'%(base, base_fraction_list_dict[base][position]*100, base_count_list_dict[base][position]) 
                    for base in NORMAL_DNA_BASES]
            display_pos = position+1 if average_both_sides else _relative_position_vs_cut(position, distance_from_middle)
            all_lines +=  " - position %s: \t%s\n"%(display_pos, ', \t'.join(data))
    return all_lines


def base_fraction_stats(base_count_position_list_dict, overall_GC_content=0.5, print_single_pvalues=False, print_summary=True, 
                        pvalue_cutoffs = [0.05, 1e-10, 1e-99], cutoff_marks=['*', '**', '***']):
    """ Given the base counts at each position, give p-values for whether they're different from the overall GC content.

    Base_count_position_list_dict should be the output of base_count_dict.

    Statistical method: according to the Handbook of Biological Statistics, what we want is a goodness-of-fit test 
      of the results vs the expected distribution (i.e. the GC content) - exact test, G-test, or Chi-square test.  
     Scipy has the chi-square test, so we're using that.  (MAYBE-TODO could also get more GoF tests from statsmodels - 
      http://statsmodels.sourceforge.net/stable/stats.html#goodness-of-fit-tests-and-measures.)

    Optionally print details and/or summary, based on the pvalue cutoffs given.
    """
    if not pvalue_cutoffs==sorted(pvalue_cutoffs, reverse=True):
        raise ValueError("pvalue_cutoffs must be sorted, largest first!")

    lengths = set([len(l) for l in base_count_position_list_dict.values()])
    if len(lengths)>1:  raise ValueError("Different bases have different count list lengths! %s"%lengths)
    length = lengths.pop()

    expected_base_fractions = base_fractions_from_GC_content(overall_GC_content)

    raw_position_pvalues = []
    FDRadj_position_pvalues = []
    base_fractions_by_pos = []
    for position,base_counts in enumerate(zip(*[base_count_position_list_dict[base] for base in NORMAL_DNA_BASES])):
        base_total = sum(base_counts)
        base_fractions = [count/base_total for count in base_counts]
        base_fractions_by_pos.append(dict(zip(NORMAL_DNA_BASES,base_fractions)))
        expected_base_fractions_list = [expected_base_fractions[base] for base in NORMAL_DNA_BASES]
        pvalue = statistics_utilities.chisquare_goodness_of_fit(base_counts, expected_base_fractions_list)
        raw_position_pvalues.append(pvalue)
    # adjust p-values for multiple testing - although it's not clear this is really needed, 
    #  since we EXPECT the significant parts to be right around the cut site, we're only checking a longer region just in case,
    #  and how long a region we're checking is pretty arbitrary...
    FDRadj_position_pvalues = statistics_utilities.FDR_adjust_pvalues(raw_position_pvalues, method='BH')

    if print_single_pvalues or print_summary:
        def base_fractions_string(base_fraction_list_dict):
            return ', '.join(['%.2f %s'%(base_fraction_list_dict[base], base) for base in NORMAL_DNA_BASES])
        print "expected base fractions: %s"%base_fractions_string(expected_base_fractions)
    if print_single_pvalues:
        relative_pos = lambda pos: _relative_position_vs_cut(pos, length/2)
        # print info for only the LOWEST cutoff matched by the pvalue
        print "single positions with raw p-value <= %s:"%max(pvalue_cutoffs)
        for position,(pvalue, adj_pvalue, base_fractions) in enumerate(zip(raw_position_pvalues, FDRadj_position_pvalues, 
                                                                           base_fractions_by_pos)):
            for cutoff,mark in reversed(zip(pvalue_cutoffs, cutoff_marks)):
                if pvalue <= cutoff:
                    print " %s pvalue %.2g (FDR-adjusted %.2g) for base %s (base fractions %s)" % (mark, pvalue, adj_pvalue, 
                                                           relative_pos(position), base_fractions_string(base_fractions))
                    break
    if print_summary:
        # grab the counts of raw and adjusted p-values over cutoffs (CUMULATIVE - a pvalue of 0 is counted for all cutoffs)
        raw_pvalue_cutoff_counts = defaultdict(lambda: 0)
        adj_pvalue_cutoff_counts = defaultdict(lambda: 0)
        for position,(pvalue, adj_pvalue) in enumerate(zip(raw_position_pvalues, FDRadj_position_pvalues)):
            for cutoff,mark in zip(pvalue_cutoffs, cutoff_marks):
                if adj_pvalue <= cutoff:     adj_pvalue_cutoff_counts[cutoff] += 1
                if pvalue <= cutoff:         raw_pvalue_cutoff_counts[cutoff] += 1
        def pvalue_cutoff_count_list(pvalue_count_dict, cutoffs):
            return ', '.join(["%s <= %s"%(pvalue_count_dict[cutoff], cutoff) for cutoff in cutoffs])
        print "out of %s positions:\n raw p-values: %s\n FDR-adjusted p-values: %s" % (length, 
                                                           pvalue_cutoff_count_list(raw_pvalue_cutoff_counts, pvalue_cutoffs), 
                                                           pvalue_cutoff_count_list(adj_pvalue_cutoff_counts, pvalue_cutoffs))
    return raw_position_pvalues, FDRadj_position_pvalues


def base_fraction_stats_compare(base_count_position_list_dict_1, base_count_position_list_dict_2, name1='SET1', name2='SET2', 
                                print_single_pvalues=False, print_summary=True, 
                                pvalue_cutoffs = [0.05, 1e-10, 1e-99], cutoff_marks=['*', '**', '***']):
    """ Like base_fraction_stats, but compare two base-count datasets to each other at east position, instead of to a GC-content. """
    if not pvalue_cutoffs==sorted(pvalue_cutoffs, reverse=True):
        raise ValueError("pvalue_cutoffs must be sorted, largest first!")

    lengths = set([len(l) for l in base_count_position_list_dict_1.values() + base_count_position_list_dict_2.values()])
    if len(lengths)>1:  raise ValueError("Different bases have different count list lengths! %s"%lengths)
    length = lengths.pop()

    raw_position_pvalues = []
    FDRadj_position_pvalues = []
    base_fractions_by_pos_1, base_fractions_by_pos_2 = [], []
    base_counts_list_1 = zip(*[base_count_position_list_dict_1[base] for base in NORMAL_DNA_BASES])
    base_counts_list_2 = zip(*[base_count_position_list_dict_2[base] for base in NORMAL_DNA_BASES])
    for base_counts_1,base_counts_2 in zip(base_counts_list_1,base_counts_list_2):
        base_total_1 = sum(base_counts_1)
        base_fractions_1 = [count/base_total_1 for count in base_counts_1]
        base_fractions_by_pos_1.append(dict(zip(NORMAL_DNA_BASES,base_fractions_1)))
        base_total_2 = sum(base_counts_2)
        base_fractions_2 = [count/base_total_2 for count in base_counts_2]
        base_fractions_by_pos_2.append(dict(zip(NORMAL_DNA_BASES,base_fractions_2)))
        pvalue = statistics_utilities.chisquare_independence(base_counts_1, base_counts_2)
        raw_position_pvalues.append(pvalue)
    # adjust p-values for multiple testing - although it's not clear this is really needed, 
    #  since we EXPECT the significant parts to be right around the cut site, we're only checking a longer region just in case,
    #  and how long a region we're checking is pretty arbitrary...
    FDRadj_position_pvalues = statistics_utilities.FDR_adjust_pvalues(raw_position_pvalues, method='BH')

    if print_single_pvalues or print_summary:
        def base_fractions_string(base_fraction_list_dict):
            return ', '.join(['%.2f %s'%(base_fraction_list_dict[base], base) for base in NORMAL_DNA_BASES])
    if print_single_pvalues:
        relative_pos = lambda pos: _relative_position_vs_cut(pos, length/2)
        # print info for only the LOWEST cutoff matched by the pvalue
        print "single positions with raw p-value <= %s:"%max(pvalue_cutoffs)
        for position,(pvalue, adj_pvalue, base_fractions_1, base_fractions_2) in enumerate(zip(raw_position_pvalues, 
                                                   FDRadj_position_pvalues, base_fractions_by_pos_1, base_fractions_by_pos_2)):
            for cutoff,mark in reversed(zip(pvalue_cutoffs, cutoff_marks)):
                if pvalue <= cutoff:
                    print " %s pvalue %.2g (FDR-adjusted %.2g) for base %s (%s %s; %s %s)" % (mark, 
                                                                                    pvalue, adj_pvalue, relative_pos(position), 
                                                                                    name1, base_fractions_string(base_fractions_1), 
                                                                                    name2, base_fractions_string(base_fractions_2))
                    break
    if print_summary:
        # grab the counts of raw and adjusted p-values over cutoffs (CUMULATIVE - a pvalue of 0 is counted for all cutoffs)
        raw_pvalue_cutoff_counts = defaultdict(lambda: 0)
        adj_pvalue_cutoff_counts = defaultdict(lambda: 0)
        for position,(pvalue, adj_pvalue) in enumerate(zip(raw_position_pvalues, FDRadj_position_pvalues)):
            for cutoff,mark in zip(pvalue_cutoffs, cutoff_marks):
                if adj_pvalue <= cutoff:     adj_pvalue_cutoff_counts[cutoff] += 1
                if pvalue <= cutoff:         raw_pvalue_cutoff_counts[cutoff] += 1
        def pvalue_cutoff_count_list(pvalue_count_dict, cutoffs):
            return ', '.join(["%s <= %s"%(pvalue_count_dict[cutoff], cutoff) for cutoff in cutoffs])
        print "out of %s positions:\n raw p-values: %s\n FDR-adjusted p-values: %s" % (length, 
                                                           pvalue_cutoff_count_list(raw_pvalue_cutoff_counts, pvalue_cutoffs), 
                                                           pvalue_cutoff_count_list(adj_pvalue_cutoff_counts, pvalue_cutoffs))
    return raw_position_pvalues, FDRadj_position_pvalues
    # TODO should probably refactor this to re-use some code from base_fraction_stats, there's a lot of duplication...


def base_fraction_plot(base_count_position_list_dict, flank_size=10, 
                       normalize_to_GC_contents=1, overall_GC_content=0.5, genome_info='', 
                       add_markers=True, bases_plotstyles={'A': 'g^-', 'T':'rv-', 'C':'bs-', 'G':'yo-'}):
    """ Plot the base fractions at each position, with given flanksize, normalized to GC content or not.

    Base_count_position_list_dict should be the output of base_count_dict.
    Normalize_to_GC_contents can be 0 (no normalization), 1 (difference between real and expected base contents), or 2 (ratio).
    """
    real_base_fraction_list_dict = base_fraction_dict_from_count_dict(base_count_position_list_dict)
    pos_after_insertion = int(len(real_base_fraction_list_dict['A']) / 2)
    expected_base_fractions = base_fractions_from_GC_content(overall_GC_content)
    for base in NORMAL_DNA_BASES:
        raw_plot_data = real_base_fraction_list_dict[base][pos_after_insertion-flank_size:pos_after_insertion+flank_size] 
        assert len(raw_plot_data) == flank_size*2
        if normalize_to_GC_contents==0:         plot_data = raw_plot_data
        elif normalize_to_GC_contents==1:       plot_data = [x-expected_base_fractions[base] for x in raw_plot_data]
        elif normalize_to_GC_contents==2:       plot_data = [x/expected_base_fractions[base] for x in raw_plot_data]
        else:                           raise Exception("normalize_to_GC_contents must be 0/1/2, not %s!"%normalize_to_GC_contents)
        if add_markers:     mplt.plot(plot_data, bases_plotstyles[base], label=base, markeredgecolor='none')
        else:               mplt.plot(plot_data, bases_plotstyles[base][0], label=base)
    mplt.legend(loc=2, prop=FontProperties(size='smaller'))
    ylabel = 'fraction of bases in given position'
    if normalize_to_GC_contents==1:     ylabel += ',\nas a difference from %s GC content'%genome_info
    elif normalize_to_GC_contents==2:   ylabel += ',\nas a ratio to %s GC content'%genome_info
    else:                               ylabel = 'raw ' + ylabel
    mplt.ylabel(ylabel,ha='center')
    # change the xticks to use -1 before the insertion position and 1 after, no 0
    xticks = range(flank_size*2)
    mplt.xlim(0,flank_size*2-1)
    mplt.xticks(xticks, [_relative_position_vs_cut(x, flank_size) for x in xticks])
    mplt.xlabel('relative genome position (dotted line is the insertion position)')
    # put a dashed line at the insertion position
    ylim = mplt.ylim()
    mplt.vlines(flank_size-0.5, *ylim, linestyles='dashed')
    mplt.ylim(*ylim)
    # MAYBE-TODO add stars/questionmarks to columns based on pvalues
    # MAYBE-TODO add info about how many sequences were in each dataset


### UNIT-TESTS

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__grab_flanking_regions(self):
        # standard inputs
        mutantfile = 'test_data/INPUT_mutants_for_flanking-regions.txt'
        test_genome = {'chr1':'AAAGGGCCC', 'chr2':'CGCG'}
        args = (mutantfile, test_genome)
        # if there's a both-strand mutant and ignore_both_strand_mutants is False, raise exception 
        self.assertRaises(ValueError, grab_flanking_regions_from_mutantfile, *args, flanksize=2, ignore_both_strand_mutants=False)
        # alwys use ignore_both_strand_mutants=True from now on, since otherwise we get an error
        kwargs = dict(ignore_both_strand_mutants=True)
        # correct flanking regions and readcounts
        assert grab_flanking_regions_from_mutantfile(*args, flanksize=0, **kwargs) == [('',10), ('',10), ('',1), ('',10)]
        assert grab_flanking_regions_from_mutantfile(*args, flanksize=1, **kwargs) == [('AG',10), ('CT',10), ('GG',1), ('GC',10)]
        assert grab_flanking_regions_from_mutantfile(*args, flanksize=2, **kwargs) == [('AAGG',10), ('CCTT',10), 
                                                                                       ('AGGG',1), ('CGCG',10)]
        # when the flanksize gets high enough that some mutants run into chromosome edges, it's padded with . 
        assert grab_flanking_regions_from_mutantfile(*args, flanksize=4, **kwargs) == [('.AAAGGGC',10), ('GCCCTTT.',10), 
                                                                                       ('AAAGGGCC',1), ('..CGCG..',10)]
        # testing that min_readcount works
        for M in (2,3,5,10):
            assert grab_flanking_regions_from_mutantfile(*args, flanksize=2, min_readcount=M, **kwargs) == [('AAGG',10), 
                                                                                                            ('CCTT',10), ('CGCG',10)]
        for M in (11, 12, 100, 1000):
            assert grab_flanking_regions_from_mutantfile(*args, flanksize=2, min_readcount=M, **kwargs) == []
        # testing that chromosome_check_function works
        f_chr1 = lambda x: x.endswith('1')
        f_chr2 = lambda x: x.endswith('2')
        f_chr3 = lambda x: x.endswith('3')
        assert grab_flanking_regions_from_mutantfile(*args,flanksize=1, chromosome_check_function=f_chr1, **kwargs) == [('AG',10),
                                                                                                                ('CT',10),('GG',1)]
        assert grab_flanking_regions_from_mutantfile(*args,flanksize=1, chromosome_check_function=f_chr2, **kwargs) == [('GC',10)]
        assert grab_flanking_regions_from_mutantfile(*args,flanksize=1, chromosome_check_function=f_chr3, **kwargs) == []

    def test__filter_flanking_regions_by_pattern(self):
        # all flanking regions must be same length
        self.assertRaises(ValueError, filter_flanking_regions_by_pattern, [('AA',2), ('AAAA',1)], '', False)
        # flanking region and pattern length must be even
        self.assertRaises(ValueError, filter_flanking_regions_by_pattern, [('AAA',2)], '', False)
        self.assertRaises(ValueError, filter_flanking_regions_by_pattern, [('',1)], 'ANN', False)
        # pattern can't be longer than flanking regions
        self.assertRaises(ValueError, filter_flanking_regions_by_pattern, [('AAA',2)], 'ANNNNN', False)
        # empty list always gives empty list
        for pattern in 'NN AAAAAA ATCG'.split():
            for orient in (True,False):
                assert filter_flanking_regions_by_pattern([], pattern, orient, False) == []
        # empty/all-N pattern matches everything
        all_4bp_regions = [(''.join(bases),1) for bases in itertools.product(*[NORMAL_DNA_BASES for _ in range(4)])]
        for pattern in ['', 'NN', 'NNNN']:
            assert filter_flanking_regions_by_pattern(all_4bp_regions, pattern, False, False) == (all_4bp_regions, [])
        # function for easier checking of more complex cases while ignoring counts
        def _check_patterns_all_counts_1(input_seqs_str, pattern, both_orient, expected_match_seqs_str):
            full_input = [(seq,1) for seq in input_seqs_str.split()]
            expected_output_match = [(seq,1) for seq in expected_match_seqs_str.split()]
            expected_match_seqs_set = set(expected_match_seqs_str.split())
            expected_nomatch_seqs = [seq for seq in input_seqs_str.split() if seq not in expected_match_seqs_set]
            expected_output_nomatch = [(seq,1) for seq in expected_nomatch_seqs]
            real_output_match, real_output_nomatch = filter_flanking_regions_by_pattern(full_input, pattern, both_orient, False)
            if not real_output_match == expected_output_match:
                print real_output_match, expected_output_match
                return False
            if not real_output_nomatch == expected_output_nomatch:
                print real_output_nomatch, expected_output_nomatch
                return False
            return True
        # a few more complex cases
        assert _check_patterns_all_counts_1('AAAT AATT GTTT', 'ANTN', False, 'AATT')
        # including a shorter pattern that needs to be padded, and seqs with Ns
        assert _check_patterns_all_counts_1('AAGT TAGN GTTT GNNT GNGT GNAT', 'AG', False, 'AAGT TAGN GNNT GNGT')
        # check that counts are propagated properly
        assert filter_flanking_regions_by_pattern([('AA',1),('AG',2),('AC',100)], 'AN', False, False)[0] == [('AA',1), 
                                                                                                             ('AG',2), ('AC',100)]
        ### Two cases for both_orientations=True!  
        ###  That one's randomized, so it's harder to test - I'm doing 100 repeats and making sure all the valid results show up.
        # 1) non-palindrome pattern:
        #  GGTT doesn't match NTTN in either direction, so it'll show up in nomatch either forward or reverse-complement.
        #  AAAT matches only when reverse-complement, and ATTG only when forward, so they'll show up as ATTT and ATTG, always.
        all_match, all_nomatch = set(), set()
        for _ in range(100):
            match,nomatch = filter_flanking_regions_by_pattern([('AAAT',1),('ATTG',2),('GGTT',3)], 'NTTN', True, False)
            all_match.add(tuple(match))
            all_nomatch.add(tuple(nomatch))
        assert all_match == set([ (('ATTT',1),('ATTG',2)) ])
        assert all_nomatch == set([ (('GGTT',3),), (('AACC',3),) ])
        # 2) palindrome pattern:    CTAC matches in either direction, so it'll show up as CTAC or GTAG;
        #                           GTTT doesn't match in either direction, so it'll show up as GTTT or AAAC.
        all_match, all_nomatch = set(), set()
        for _ in range(100):
            match,nomatch = filter_flanking_regions_by_pattern([('CTAC',1),('GTTT',2)], 'NTAN', True, False)
            all_match.add(tuple(match))
            all_nomatch.add(tuple(nomatch))
        assert all_match   == set([ (('CTAC',1),), (('GTAG',1),) ])
        assert all_nomatch == set([ (('GTTT',2),), (('AAAC',2),) ])

    def test___relative_position_vs_cut(self):
        assert _relative_position_vs_cut(1,2) == -1
        assert _relative_position_vs_cut(2,2) == 1
        for N in range(1000):
            assert _relative_position_vs_cut(N,0) == N+1
            assert _relative_position_vs_cut(N,N) == 1
            assert _relative_position_vs_cut(N-1,N) == -1

    def test__print_base_count_fraction_for_dist(self):
        # these tests are pretty brief, since the function returns strings and is annoying to test
        # fails if lengths aren't all the same, or seq or ignore_bases_pattern aren't even, 
        #  or if ignore_bases_pattern is longer than 2*distance_from_middle
        self.assertRaises(ValueError, print_base_count_fraction_for_dist, [('AT',1), ('GGG',3)], 1)
        self.assertRaises(ValueError, print_base_count_fraction_for_dist, [('ATG',1)], 1)
        self.assertRaises(ValueError, print_base_count_fraction_for_dist, [('AT',1)], 1, ignore_bases_pattern='T')
        self.assertRaises(ValueError, print_base_count_fraction_for_dist, [('AT',1)], 1, ignore_bases_pattern='TTTT')
        # basic functionality
        assert print_base_count_fraction_for_dist([('AT',1), ('GG',3)], 1, ignore_bases_pattern=None, average_both_sides=False) == ( 
            ' - position -1: \tA 25% (1), \tC 0% (0), \tT 0% (0), \tG 75% (3)\n'
           +' - position 1: \tA 0% (0), \tC 0% (0), \tT 25% (1), \tG 75% (3)\n')
        # ignore_bases_pattern
        assert print_base_count_fraction_for_dist([('AT',1), ('GG',3)], 1, ignore_bases_pattern='AA', average_both_sides=False) == ''
        assert print_base_count_fraction_for_dist([('AT',1), ('GG',3)], 1, ignore_bases_pattern='AN', average_both_sides=False) == ( 
            ' - position 1: \tA 0% (0), \tC 0% (0), \tT 25% (1), \tG 75% (3)\n')
        # average_both_sides
        assert print_base_count_fraction_for_dist([('AT',1), ('GG',1)], 1, ignore_bases_pattern=None, average_both_sides=True) == ( 
            ' - position 1: \tA 0% (0), \tC 25% (1), \tT 50% (2), \tG 25% (1)\n')
        assert print_base_count_fraction_for_dist([('AT',1), ('GG',3)], 1, ignore_bases_pattern=None, average_both_sides=True) == ( 
            ' - position 1: \tA 0% (0), \tC 38% (3), \tT 25% (2), \tG 38% (3)\n')
        # average_both_sides AND ignore_bases_pattern, since they interact in complicated ways!
        assert print_base_count_fraction_for_dist([('AT',1), ('GG',3)], 1, ignore_bases_pattern='AA', average_both_sides=True) == ''
        assert print_base_count_fraction_for_dist([('CATG',1)], 2, ignore_bases_pattern='ANNA', average_both_sides=True) == ( 
            ' - position 1: \tA 0% (0), \tC 0% (0), \tT 100% (2), \tG 0% (0)\n')
        assert print_base_count_fraction_for_dist([('CATG',1)], 2, ignore_bases_pattern='NAAN', average_both_sides=True) == ( 
            ' - position 2: \tA 0% (0), \tC 0% (0), \tT 0% (0), \tG 100% (2)\n')
        assert print_base_count_fraction_for_dist([('CATG',1)], 2, ignore_bases_pattern='ANAN', average_both_sides=True) == ( 
            ' - position 1: \tA 0% (0), \tC 0% (0), \tT 100% (1), \tG 0% (0)\n'
           +' - position 2: \tA 0% (0), \tC 0% (0), \tT 0% (0), \tG 100% (1)\n')

    # LATER-TODO add unit-tests for all the other functions?  Lower priority, since they're either straightforward, or mostly focus on printing/plotting data rather than complicated transformations.

if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
