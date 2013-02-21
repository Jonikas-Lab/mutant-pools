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
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
R_stats = importr('stats')
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties
# my modules
import mutant_analysis_classes
from basic_seq_utilities import get_all_seq_length, base_count_dict, base_fraction_dict, base_fraction_dict_from_count_dict, base_fractions_from_GC_content, NORMAL_DNA_BASES, reverse_complement, check_seq_against_pattern, write_fasta_line
import general_utilities


def grab_flanking_regions(mutant_dataset_infile, genome, flanksize=200, 
                          min_readcount=0, chromosome_check_function=lambda x: True, ignore_both_strand_mutants=False):
    """ Return (flanking_seq,readcount) with both-side genomic flanking sequences for insertional mutants in mutant_dataset_infile.

    Grab all the insertion positions from mutant_dataset_infile (pickled mutant_analysis_classes.Insertional_mutant_dataset object), 
     use genome (a chrom_name:seq dict) to figure out the flanksize-length flanking sequences on both sides
      (padded with '.' if the end of the chromosome is too close), and write them to outfile in fasta format. 
    The flanking sequences will be in the same orientation as the insertion cassette.

    Filter the mutants: 
     - by readcount - ignore mutants with total readcount below min_readcount=0
     - by chromosome - ignore mutants in chromosomes for which chromosome_check_function returns False
     - by strand - both-strand (merged tandem) mutants will be ignored if ignore_both_strand_mutants is True, 
                otherwise ValueError will be raised; ValueError will be raised for other unexpected strand values.
    For all remaining mutants, append (flanking region seq, total_readcount) to output list.
    """
    flanking_region_count_list = []
    count_all, count_lowreads, count_midreads, count_highreads = 0, 0, 0, 0
    dataset = general_utilities.unpickle(mutant_dataset_infile)
    for mutant in dataset:
        # filter out mutants with wrong readcounts or in wrong chromosomes
        if not chromosome_check_function(mutant.position.chromosome):           continue
        if mutant.total_read_count < min_readcount:                             continue
        # filter out both-stranded mutants if desired; 
        if mutant.position.strand not in '+-':  
            if mutant.position.strand=='both' and ignore_both_strand_mutants:   continue
            else:                                   raise ValueError("Unexpected mutant strand! %s"%mutant.position)
        # grab mutant position/chromosome
        pos_before_insertion = mutant.position.min_position
        chromosome_seq = genome[mutant.position.chromosome]
        # ignore cassette tandems (i.e. insertions that map to start or end of cassette)
        if mutant_analysis_classes.is_cassette_chromosome(mutant.position.chromosome):
            if pos_before_insertion in [0, len(chromosome_seq)]:  continue
        # grab positions to get the sequence flanksize bp on both sides of the insertion;
        #  pad it with . if there's a chromosome start/end closer than that.
        flank_start = pos_before_insertion - flanksize
        flank_end = pos_before_insertion + flanksize
        start_padding = '.' * max(0, -flank_start)
        flank_start = max(flank_start, 0)
        end_padding = '.' * max(0, flank_end - len(chromosome_seq))
        flank_end = min(flank_end, len(chromosome_seq))
        # grab the actual flanking sequence, add padding; 
        #  reverse-complement it if mutant is -strand, to keep them in the same orientation as the cassette insertion
        full_flanking_seq = chromosome_seq[flank_start:flank_end]
        full_flanking_seq = start_padding + full_flanking_seq + end_padding
        if mutant.position.strand == '-':
            full_flanking_seq = reverse_complement(full_flanking_seq)
        # append the sequence and readcount to output data
        flanking_region_count_list.append((full_flanking_seq, mutant.total_read_count))
    return flanking_region_count_list


def filter_flanking_regions_by_pattern(flanking_region_count_list, pattern, either_orientation=True):
    """ Return separate lists of flanking regions that do and don't match given sequence pattern.
    
    flanking_region_count_list should be a list of (flanking_region, count) pairs (like from grab_flanking_regions); 
     the two return values (flanking regions that match and don't match the pattern) are the same format.

    The pattern should be a sequence string (allowed letters are ACTGN). It'ss considered to be centered around the cut site; 
     the flanking regions likewise.  E.g. if pattern is GNAN, a fl.region of GCAC or TTGCACTT would match, but TTTTGCAC would not. 
    If either_orientation is True, each flanking region will be tried against the pattern in both the forward and the reverse
     orientation, and the returned flanking region will be in the orientation that matched - e.g. if pattern is GNAN, 
     a flanking region of either TTGCACTT or TTCTCCTT would match (forward and rev-compl respectively), 
      and the latter would be returned as rev-compl, AAGGAGAA.
    """
    if not flanking_region_count_list:  return []
    flanking_region_length = get_all_seq_length(zip(*flanking_region_count_list)[0])
    if flanking_region_length % 2:              sys.exit("Flanking region length must be an even number!")
    if len(pattern) % 2:                        sys.exit("Pattern length must be an even number!")
    if len(pattern) > flanking_region_length:   sys.exit("Pattern cannot be longer than flanking regions!")
    # pad the pattern to match the flanking region length
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
    return flanking_region_count_list_match, flanking_region_count_list_nomatch
    # TODO amend unit-tests to reflect the randomization and other changes!


def print_flanking_regions_to_fasta(flanking_region_count_list, outfile, convert_counts=lambda x: x):
    """ Given a (seq,count) list, make fasta file with the each seq present convert_counts(count) times. """
    with open(outfile, 'w') as OUTFILE:
        for N,(seq,count) in enumerate(flanking_region_count_list):
            seqname = '%s (%s reads)' % (N, count)
            for _ in range(convert_counts(count)):
                write_fasta_line(seqname, seq, OUTFILE)


def _relative_position_vs_cut(pos, reference_pos):
    """ Convert zero-based position into relative position around a cut site after reference_pos, with no 0.

    If refpos is 1, we assume the true reference cut position is between 1 and 2, so positions 1 and 2 become -1 and 1, 
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
    if flanking_region_length % 2:              sys.exit("Flanking region length must be an even number!")
    flank_length = int(flanking_region_length/2)
    # grab only the flanking region length we're actually interested in, and convert the counts
    local_flanking_region_length = 2*distance_from_middle
    local_flanking_region_count_list = [(seq[flank_length-distance_from_middle:flank_length+distance_from_middle], 
                                         convert_counts(count)) for (seq,count) in flanking_region_count_list]
    # apply ignore_bases_pattern by changing the relevant bases to N
    if ignore_bases_pattern is not None:
        if len(ignore_bases_pattern) % 2:       sys.exit("Ignore_bases_pattern length must be an even number!")
        length_diff = int((len(ignore_bases_pattern)-local_flanking_region_length)/2)
        if length_diff>0:   ignore_bases_pattern = ignore_bases_pattern[length_diff:-length_diff]
        if length_diff<0:   ignore_bases_pattern = 'N'*length_diff + ignore_bases_pattern + 'N'*length_diff
        def mask_seq(seq, mask_pattern=ignore_bases_pattern):
            return ''.join([(N if if_mask=='N' else base) for (base,if_mask) in zip(seq, mask_pattern)])
        local_flanking_region_count_list = [(mask_seq(seq), count) for  (seq,count) in local_flanking_region_count_list]
    # if average_both_sides, make a new local_flanking_region_count_list that has each half of each sequence separately
    if average_both_sides:
        new_flanking_region_count_list = []
        for flanking_region,count in local_flanking_region_count_list:
            first_half = reverse_complement(flanking_region[:distance_from_middle])
            second_half = flanking_region[distance_from_middle:]
            new_flanking_region_count_list.extend([(first_half,count), (second_half,count)])
        local_flanking_region_count_list = new_flanking_region_count_list
    base_count_list_dict = base_count_dict(local_flanking_region_count_list)
    base_fraction_list_dict = base_fraction_dict(local_flanking_region_count_list)
    all_lines = ''
    for position in range(local_flanking_region_length):
        display_pos = position if average_both_sides else _relative_position_vs_cut(position, distance_from_middle)
        data = ['%s %.0f%% (%s)'%(base, base_fraction_list_dict[base][position]*100, base_count_list_dict[base][position]) 
                for base in NORMAL_DNA_BASES]
        all_lines +=  " - position %s: \t%s\n"%(display_pos, ', \t'.join(data))
    return all_lines
    # TODO unit-test! This got rather complicated...


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
        sys.exit("pvalue_cutoffs must be sorted, largest first!")

    lengths = set([len(l) for l in base_count_position_list_dict.values()])
    if len(lengths)>1:  sys.exit("Different bases have different count list lengths! %s"%lengths)
    length = lengths.pop()

    expected_base_fractions = base_fractions_from_GC_content(overall_GC_content)

    raw_position_pvalues = []
    FDRadj_position_pvalues = []
    base_fractions_by_pos = []
    for position,base_counts in enumerate(zip(*[base_count_position_list_dict[base] for base in NORMAL_DNA_BASES])):
        base_total = sum(base_counts)
        base_fractions = [count/base_total for count in base_counts]
        base_fractions_by_pos.append(dict(zip(NORMAL_DNA_BASES,base_fractions)))
        # NOTE - for scipy.stats.chisquare the reference has to be normalized so it adds up to the same total! WEIRD. So do that. 
        expected_base_counts = [expected_base_fractions[base]*base_total for base in NORMAL_DNA_BASES]
        # Inputs have to be as scipy.array.  The return value is a (chisquare_statistic, pvalue) tuple - just grab the pvalue.
        pvalue = scipy.stats.chisquare(scipy.array(base_counts), scipy.array(expected_base_counts))[1]
        raw_position_pvalues.append(pvalue)
    # adjust p-values for multiple testing - although it's not clear this is really needed, 
    #  since we EXPECT the significant parts to be right around the cut site, we're only checking a longer region just in case,
    #  and how long a region we're checking is pretty arbitrary...
    FDRadj_position_pvalues = R_stats.p_adjust(FloatVector(raw_position_pvalues), method='BH')

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


def base_fraction_plot(base_count_position_list_dict, flank_size=10, 
                       normalize_to_GC_contents=1, overall_GC_content=0.5, 
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
    if normalize_to_GC_contents==1:     ylabel += ',\nas a difference from overall GC content of that genome'
    elif normalize_to_GC_contents==2:   ylabel += ',\nas a ratio to overall GC content of that genome'
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

    def test__filter_flanking_regions_by_pattern(self):
        # all flanking regions must be same length
        self.assertRaises(SystemExit, filter_flanking_regions_by_pattern, [('AA',2), ('AAAA',1)], '')
        # flanking region and pattern length must be even
        self.assertRaises(SystemExit, filter_flanking_regions_by_pattern, [('AAA',2)], '')
        self.assertRaises(SystemExit, filter_flanking_regions_by_pattern, [('',1)], 'ANN')
        # pattern can't be longer than flanking regions
        self.assertRaises(SystemExit, filter_flanking_regions_by_pattern, [('AAA',2)], 'ANNNNN')
        # empty list always gives empty list
        for pattern in 'NN AAAAAA ATCG'.split():
            for orient in (True,False):
                assert filter_flanking_regions_by_pattern([], pattern, either_orientation=orient) == []
        # empty/all-N pattern matches everything
        all_4bp_regions = [(''.join(bases),1) for bases in itertools.product(*[NORMAL_DNA_BASES for _ in range(4)])]
        for pattern in ['', 'NN', 'NNNN']:
          for orient in (True,False):
            assert filter_flanking_regions_by_pattern(all_4bp_regions, pattern, orient) == all_4bp_regions
        assert filter_flanking_regions_by_pattern
        # function for easier checking of more complex cases while ignoring counts
        def _check_patterns_all_counts_1(input_seqs_str, pattern, both_orient, output_seqs_str):
            full_input = [(seq,1) for seq in input_seqs_str.split()]
            expected_output = [(seq,1) for seq in output_seqs_str.split()]
            real_output = filter_flanking_regions_by_pattern(full_input, pattern, both_orient)
            return (real_output == expected_output)
        # a few more complex cases
        assert _check_patterns_all_counts_1('AAAT AATT GTTT', 'ANTN', False, 'AATT')
        assert _check_patterns_all_counts_1('AAAT AATT GTTT', 'ANTN', True, 'ATTT AATT')
        # including a shorter pattern that needs to be padded, and seqs with Ns
        assert _check_patterns_all_counts_1('AAGT TAGN GTTT GNNT GNGT GNAT', 'AG', False, 'AAGT TAGN GNNT GNGT')
        # TODO finish!

        # check that counts are propagated properly
        assert filter_flanking_regions_by_pattern([('AA',1), ('AG',2), ('AC',100)], 'AN', False) == [('AA',1), ('AG',2), ('AC',100)]

    # TODO write unit-tests for some of the other more complicated functions!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
