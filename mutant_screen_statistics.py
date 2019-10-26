#! /usr/bin/env python2.7

"""
Basic analysis pipeline for pooled mutant screens, resulting in p-values for each gene based on allele numbers and phenotypes: 
    - filter IBs by minimum control readcount; calculate sample:control phenotype ratios for each IB, 
    - sort mutants into bins based on the ratios and cutoffs (non-hits, weak/medium/strong hits)
    - for each gene, do a Fisher's exact test comparing the numbers of mutants in this gene in each bin to all mutants, 
        optionally filtering by gene feature and confidence level
    - do FDR-correction on the resulting p-values, optionally excluding genes with <N alleles
    - print results to text file, and generate a python pickle file.
    - generate scatterplot of raw and normalized IB readcounts with given readcount and phenotype cutoffs.

Additionally, if two sample names are given (usually makes sense for replicates):
    - do the above analysis for each sample, printing a single text file
    - generate a scatterplot of gene FDR values between the two samples.
If two control names are given, each control is paired with the corresponding sample for analysis; if a single control is given with two samples, both samples are analyzed with the same control.

Input file should be a python pickle file containing a {sample_name:{IB:(raw_readcount, normalized_readcount)}} dictionary, and all sample and control names provided in the options must be present in the dictionary.

The code is largely based on a simpler analysis of earlier screen data done for Xiaobo in ../../../arrayed_library/1603_large-lib_figures-1/notes.txt section "### do statistics to find potential gene hits"; the more complex analysis method (with multiple thresholds for Fisher's exact test) is based on the analysis Robert did for the large-scale screen paper (2018?) - see code in ../../1802_Moshe+Frieder_new-screens/Robert_pipeline, although his code isn't very readable for me and I didn't use it.

 -- Weronika Patena, 2018

USAGE: mutant_screen_statistics.py [options] input_file output_folder
"""

# standard library
from __future__ import division
import sys
import os
import unittest
import collections
# other packages
import scipy.stats
import numpy
import matplotlib.pyplot as mplt
# my modules
import general_utilities
import mutant_utilities
import mutant_IB_RISCC_classes
import statistics_utilities
import plotting_utilities
import parse_annotation_file

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    ### functionality options
    parser.add_option('-S', '--sample', type='string', default='', metavar='"Sample1"', 
              help="Picklefile of the sample {IB:(raw,norm)} readcount dictionary, or name of the sample in the -D dictionary"
                  +" (required) - put quotes around it if it has spaces or weird characters."
                  +" Can be two comma-separated values - in that case both will be analyzed separately and then as replicates.")
    parser.add_option('-C', '--control', type='string', default='', metavar='"Control1"', 
              help="Picklefile of the control {IB:(raw,norm)} readcount dictionary, or name of the control in -D dictionary"
                  +" (required) - put quotes around it if it has spaces or weird characters."
                  +" If two --sample values were given, this can be two comma-separated values or a single value.")
    parser.add_option('-D', '--data_dictionary', type='string', default='', metavar='Data', 
              help="Python pickle file containing a {sample_name:{IB:(raw_readcount, normalized_readcount)}} dictionary,"
                  +"containing all --sample and --control names and optionally any other data."
                  +" Required only if --sample/--control values are names rather than separate picklefiles.")
    parser.add_option('-p', '--phenotype_thresholds', type='string', default='0.0625,0.125,0.25,0.5', metavar='X,Y', 
              help="List of sample/control ratio thresholds for bins, comma-separated, no spaces, lowest first (default %default)."
                  +"The real thresholds used will be made two-sided, i.e. '10' -> '0.1,10' etc, unless -O is used."
                  +"If -P, this should be a list of offsets for Poisson curves instead (integers 1-6, positive or negative).")
    parser.add_option('-P', '--Poisson_thresholds', action='store_true', default=False,
              help="Instead of linear thresholds, use Poisson curves with offsets (default: %default)")
    parser.add_option('-O', '--one_sided_thresholds', action='store_true', default=False,
              help="Only look for slower-growing and not faster-growing mutants (default: two-sided)")
    parser.add_option('-m', '--min_reads', type='int', default='50', metavar='M', 
              help="Minimum raw control reads (or control+sample reads if -M) to include IB in analysis (default %default).")
    parser.add_option('-M', '--symmetrical_min_reads', action='store_true', default=False,
              help="Minimum readcount uses control+sample reads rather than only control reads (default %default).")
    parser.add_option('-x', '--min_reads_2', type='int', default='0', metavar='M', 
              help="Minimum raw SAMPLE reads to include IB in analysis - rarely needed! (default %default).")
    parser.add_option('-c', '--max_conf', type='int', default=4, metavar='C',
              help="Only mutants with mapping confidence from 1 to C will be included in the analysis (default %default)")
    parser.add_option('-f', '--features', type='string', default="CDS,intron,5'UTR,5'UTR_intron", metavar='F,F',
              help="Only mutants in these gene features will be included in the analysis (default %default)")
    parser.add_option('-a', '--min_alleles_for_FDR', type='int', default=1, metavar='A',
              help="Only genes with at least A qualified alleles are included in FDR correction (default %default)")
    parser.add_option('-F', '--full_library', action='store_true', default=False,
              help="The screen was done on the full library (default: rearray)")
    parser.add_option('-W', '--workspace_size', type='float', default=2e5,
              help="Memory used for statistics - try increasing to 2e6 if you're getting \"LDSTP is too small for this problem.  "
                  +"Try increasing the size of the workspace\" errors, but 2e6 seems to be ~20G... (default %default)")
    return parser


########################################### Basic pipeline #################################################

def modify_phenotype_thresholds(phenotype_thresholds, Poisson=False, one_sided_thresholds=False):
    """ Make sure thresholds are in the right order; make them two-sided if needed.
    """
    phenotype_thresholds = [float(x) for x in phenotype_thresholds.split(',')]
    phenotype_thresholds = [int(x) if round(x)==x else x for x in phenotype_thresholds]
    middle = 0 if Poisson else 1
    if min(phenotype_thresholds) <= middle <= max(phenotype_thresholds):
        raise Exception("All thresholds should be >%s or all <%s! (current values are %s). "%(middle, middle, phenotype_thresholds)
                        +"Symmetrical ones will be added automatically to make the test two-sided if one_sided_thresholds is False.")
    # in one-sided cases, always make it so the wt bin is at the end - makes display better and everything else easier.
    low = -10 if Poisson else 0
    high = 10 if Poisson else float('inf')
    opposite = (lambda x: -x) if Poisson else (lambda x: 1/x)
    if one_sided_thresholds:
        if numpy.mean(phenotype_thresholds) < middle:
            phenotype_thresholds = [low] + sorted(phenotype_thresholds) + [high]
        else:
            phenotype_thresholds = [high] + sorted(phenotype_thresholds, reverse=True) + [low]
    # in two-sided cases, make it so slow-growing is at the start, for consistency/readability
    #   (in normal thresholds slow-growing is fractions; in Poisson it's positive values, so reverse the sort)
    else:
        phenotype_thresholds = sorted([low] + [opposite(x) for x in phenotype_thresholds] + phenotype_thresholds + [high], 
                                     reverse=Poisson)
    print "phenotype thresholds%s: "%(' (Poisson)' if Poisson else ''), phenotype_thresholds
    return phenotype_thresholds


def filter_IBs_by_min_readcount(screen_data, min_readcount, symmetrical_min=False, min_readcount_2=0):
    """ Modifies screen_data to remove any IBs below the min readcounts.

    Returns separate dictionary of data removed from screen_data.
    """
    if symmetrical_min and min_readcount_2:
        raise Exception("Can't have a minimum readcount on control+sample AND on sample!")
    removed_screen_data = {}
    for (IB,data) in screen_data.items():
        if not symmetrical_min:
            if data[2] < min_readcount or data[0] < min_readcount_2:
                del screen_data[IB]
                removed_screen_data[IB] = data
        else:
            if data[2]+data[0] < min_readcount:
                del screen_data[IB]
                removed_screen_data[IB] = data
    return removed_screen_data


def poisson_function(offset):
    """ Given an offset (semi-arbitrary number >0; 7=infinity), return Poisson curve function.
    """
    # this offset thing in two places is a bit weird, I just hand-picked it to look good with (1,2,3,4). 
    #   the -.2 is just so the line doesn't overlap the actual data points, to make it clearer which side they count on.
    if offset > 6: 
        return lambda x: float('inf')
    return lambda x: scipy.stats.poisson(x/1.8**offset).interval(1-1/10**(4+2*offset))[0]-.2


def check_threshold_Poisson(x, y, threshold):
    """ Check whether a given (x,y) position is above a custom Poisson curve with the given offset.
    """
    # deal with the infinity special cases separately - nothing should ever be below the -inf curve or above the +inf one:
    if threshold > 6:       return True
    elif threshold < -6:    return False
    # otherwise do a normal check for a positive number, and a reverse check (swapping x and y) for negative numbers
    elif threshold > 0:
        poisson_func = poisson_function(threshold)
        return (y > poisson_func(x))
    else:
        poisson_func = poisson_function(-threshold)
        return (x < poisson_func(y))


def bin_IBs_by_phenotype(screen_data, phenotype_thresholds, Poisson=False):
    """ Bins IBs by sample/control ratio or Poisson curve based on phenotype_thresholds; returns list of IB sets for each bin.
    """
    # screen_data is an IB:(raw_sample, norm_sample, raw_ctrl, norm_ctrl) readcount dict.
    # if normal thresholds, use the normalized data
    if not Poisson:
        # need to use a high value to deal with the x/0 case; can't use inf because then it doesn't get binned because inf !< inf.
        highval = 2*max(phenotype_thresholds[1:-1])
        if_between = lambda low, high, data: (min(low,high) <= (data[1]/data[3] if data[3] else highval) < max(low,high))
    # if Poisson, have to use the raw data; control is x axis and sample is y axis. 
    # this is SLOW, ~30min for a normal dataset, but I'm not sure I can speed up Poisson math...
    else:
        if_between = lambda low, high, data: (check_threshold_Poisson(data[2],data[0],low)
                                              and not check_threshold_Poisson(data[2],data[0],high))
    binned_IBs = [set(IB for IB,data in screen_data.items() if if_between(a,b,data))
                  for (a,b) in zip(phenotype_thresholds, phenotype_thresholds[1:])]
    print "numbers of IBs in each phenotype bin: ", [len(x) for x in binned_IBs]
    if len(set.union(*binned_IBs)) != sum(len(x) for x in binned_IBs):
        print "PROBLEM: bins are double-counting IBs?? %s IBs, bin total %s"%(len(set.union(*binned_IBs)), 
                                                                              sum(len(x) for x in binned_IBs))
    elif len(screen_data) != sum(len(x) for x in binned_IBs):
        missing_IBs = set(screen_data.keys()) - set.union(*binned_IBs)
        print "PROBLEM: %s IBs couldn't be binned! %s"%(len(missing_IBs), [(IB, screen_data[IB]) for IB in missing_IBs])
    return binned_IBs
    # there is a minor asymmetry here because the bins are defined as low <= x < high, so values exactly on the edge 
    #   will belong to one bin or the other depending on the order of the thresholds.  
    # If I really wanted to be extra-careful I could make it symmetrical by making everything go toward the bin that's closer to 
    #   1 (the middle) instead of closer to 0, but that would make the code more complicated and seems unnecessary.


def raw_screen_data_per_gene(library_data_by_IB, screen_data, features, max_conf):
    """ Make a list of screen per-IB data lines for each gene, filtering by feature and confidence.

    Return a gene:IB_data dict, where IB_data is a list of (IB, (raw_sample, norm_sample, raw_ctrl, norm_ctrl)) tuples
    """
    screen_lines_per_gene = collections.defaultdict(list)
    for IB in screen_data.keys():
        x = library_data_by_IB[IB]
        if x.gene not in mutant_IB_RISCC_classes.SPECIAL_GENE_CODES.all_codes and x.feature != 'intergenic':
            if x.confidence_level <= max_conf:
                # remember some positions are in multiple genes joined by &, same for features
                for (gene,feature) in zip(x.gene.split(' & '), x.feature.split(' & ')):
                    # if_right_feature deals with the /- and |-split features
                    if features is None or mutant_utilities.if_right_feature(feature, features):
                        gene = mutant_utilities.strip_version(gene)
                        screen_lines_per_gene[mutant_utilities.strip_version(gene)].append((IB,screen_data[IB]))
    return screen_lines_per_gene


def filter_IBs_per_gene(screen_lines_per_gene, symmetric, library_data_by_IB, if_full_library=False):
    """ Filter lists of screen data for each gene to remove multiple ones in the same mutant
    
    If symmetric is True, for each mutant keep the IB with the highest control+sample readcount; otherwise just highest control.
    screen_lines_per_gene is the output of raw_screen_data_per_gene.
    Returns gene:IB_set dictionary.
    """
    if if_full_library:     get_mutant_ID = lambda x: (x.plate, x.well)
    else:                   get_mutant_ID = lambda x: x.mutant_ID
    IBs_per_gene_filtered = {}
    for gene,screen_lines in screen_lines_per_gene.items():
        screen_lines_by_mutant = collections.defaultdict(list)
        for (IB,data) in screen_lines:
            screen_lines_by_mutant[get_mutant_ID(library_data_by_IB[IB])].append((IB,data))
        filtered_IBs = set()
        for line_set in screen_lines_by_mutant.values():
            # sort by IB after the readcounts just to make sure everything is deterministic even if two IBs have the same readcount.
            if symmetric:   line_set.sort(key = lambda x: (x[1][0]+x[1][2], x[0]))
            else:           line_set.sort(key = lambda x: (x[1][2], x[0]))
            filtered_IBs.add(line_set[-1][0])
        IBs_per_gene_filtered[gene] = filtered_IBs
    return IBs_per_gene_filtered


def get_gene_bin_counts(IBs_per_gene_filtered, binned_IBs_by_phenotype):
    """ Get counts of alleles in each hit/nonhit bin for each gene

    Genes with no alleles in any bins aren't included.
    """
    gene_bin_counts = {}
    for (gene, filtered_IBs) in IBs_per_gene_filtered.items():
        bin_counts = [len(filtered_IBs & bin_IBs) for bin_IBs in binned_IBs_by_phenotype]
        if sum(bin_counts):
            gene_bin_counts[gene] = bin_counts
    return gene_bin_counts


def gene_statistics(gene_bin_counts, min_alleles_for_FDR=1, workspace_size=2e5):
    """ Calculate a p-value (using Fisher's exact test) and FDR (BH method) for each gene compared to all alleles

    Output: gene:[binned_allele_counts, pval, FDR] dictionary, where FDR is NaN for genes with <min_alleles_for_FDR alleles.
    """
    bin_totals = [sum(x) for x in zip(*gene_bin_counts.values())]
    print "numbers of filtered IBs in each phenotype bin: ", bin_totals
    # if there are only 2 bins, can use scipy.stats.fisher_exact; otherwise use my custom one that goes through R
    #  the one using R will sometimes throw an "Increase workspace or consider using 'simulate.p.value=TRUE'" error - 
    #   for that reason I changed that to optionally use a higher workspace value, we'll see if that works.
    # MAYBE-TODO could also switch to chi-squared test for bigger numbers, but that's not great when some numbers are small...
    if len(bin_totals) == 1:    fisher_exact = lambda x: scipy.stats.fisher_exact(x)[1]
    else:                       fisher_exact = lambda x: statistics_utilities.fisher_exact(x, workspace=workspace_size)
    gene_pvals = {g:fisher_exact([bin_counts,bin_totals]) for (g,bin_counts) in gene_bin_counts.items()}
    if any(numpy.isnan(x) for x in gene_pvals.values()): raise Exception("Some pvals are NaN! Look into this!")
    # the FDR-correction has to be done on a list of pvalues, so separate out the genes that meet min_alleles_for_FDR
    genes_with_enough_alleles = sorted(g for (g,bin_counts) in gene_bin_counts.items() if sum(bin_counts) >= min_alleles_for_FDR)
    FDRs_for_some_genes = statistics_utilities.FDR_adjust_pvalues([gene_pvals[g] for g in genes_with_enough_alleles], method='BH')
    FDR_dict = collections.defaultdict(lambda: float('NaN'), zip(genes_with_enough_alleles, FDRs_for_some_genes))
    # join the bin counts, pvals and FDRs into a single dictionary, with NaN FDRs as needed
    gene_stats_data = {g: [bin_counts, gene_pvals[g], FDR_dict[g]] for (g,bin_counts) in gene_bin_counts.items()}
    return gene_stats_data


def print_mutant_read_numbers(screen_data, sample_name, control_name, header='Samples'):
    sumF = lambda x: general_utilities.format_number(sum(x), 0, 1)
    lenF = lambda x: general_utilities.format_number(len(x), 0, 1)
    overlap_IBs = set(IB for IB,x in screen_data.items() if x[0] and x[2])
    print "%s: %s %s mutants (%s reads);  Control %s %s (%s);  Overlap %s (%s/%s)"%(header,
       sample_name,  sumF(1 for x in screen_data.values() if x[0]), sumF(x[0] for x in screen_data.values()),
       control_name, sumF(1 for x in screen_data.values() if x[2]), sumF(x[2] for x in screen_data.values()),
       lenF(overlap_IBs), sumF(screen_data[x][0]  for x in overlap_IBs), sumF(screen_data[x][2] for x in overlap_IBs))


def print_ratio_counts(gene_bin_counts, wt_bins = .95, one_sided_thresholds=False):
    """ Print how many genes have wt:phenotype allele ratios of 1:0, 1:1, etc.
    
    Simplify the bins into just two bins (wt-like and phenotype, with both faster-growing and slower-growing counted as phenotype):
        if wt_bins is an integer, use the first X bins as wt (if one_sided_thresholds is True),
            or X*2-1 middle bins as wt (if one_sided_thresholds is False); 
        if it's a float<1, use enough many middle bins as wt to have wt be at least X fraction of all counts. 
    """
    if (wt_bins >= 1 and (wt_bins-1) % 2) or wt_bins < 0.3:
        raise Exception("wt_bins must be either an odd integer (number of central bins) "
                        +"or a float 0.3-0.99 (pick central bins to cover at least that fraction of mutants)")
    # first, if thresholds were two-sided (wt in the middle), "fold" them in half to get a one-sided version
    N_bins = len(gene_bin_counts.values()[0])
    middle_bin = int((N_bins-1)/2)
    if not one_sided_thresholds:
        gene_bin_counts_1side = {}
        for g,b in gene_bin_counts.items():
            b_1side = [b[i]+b[-i-1] for i in range(middle_bin)] + [b[middle_bin]] 
            assert sum(b_1side) == sum(b)
            gene_bin_counts_1side[g] = b_1side
        gene_bin_counts = gene_bin_counts_1side
    # if wt_bins defined by fraction of mutant, figure out how many are needed
    if wt_bins < 1:
        bin_totals = [sum(x) for x in zip(*gene_bin_counts.values())]
        total = sum(bin_totals)
        for tmp in range(1,len(bin_totals)):
            if sum(bin_totals[-tmp:])/total >= wt_bins:
                wt_bins = tmp
                break
        else: 
            print "(Couldn't find wt-like bins covering %s of the mutants; instead using all-but-last bins, covering %.2f.)"%(
                wt_bins, sum(bin_totals[1:])/total)
            wt_bins = N_bins-1
    # convert all-bin counts to just wt,phenotype bin counts for each gene (as a list, since gene names aren't needed)
    bin_counts_simplified = [(sum(bin_counts[:wt_bins]), sum(bin_counts[wt_bins:])) for bin_counts in gene_bin_counts.values()]
    # Count genes in different categories - should all be mutually exclusive and cover everything
    inf = float('inf')
    print_data = []
    total_between_all_cases = 0 
    for (name, min_ph, max_ph, min_wt, max_wt) in [('3+:0',   3, inf,   0, 0),
                                                   ('2:0',    2, 2,     0, 0),   
                                                   ('1:0',    1, 1,     0, 0),   
                                                   ('3+:1',   3, inf,   1, 1),   
                                                   ('2:1',    2, 2,     1, 1),   
                                                   ('1:1',    1, 1,     1, 1),   
                                                   ('2+:2+',  2, inf,   2, inf), 
                                                   ('1:2+',   1, 1,     2, inf), 
                                                   ('0:1',    0, 0,     1, 1),   
                                                   ('0:2+',   0, 0,     2, inf)]: 
        N = sum(1 for ph,wt in bin_counts_simplified if min_wt<=wt<=max_wt and min_ph<=ph<=max_ph)
        print_data.append("%s %s"%(name, N))
        total_between_all_cases += N
    print "Genes by phen:wt allele counts: " + ', '.join(print_data)
    # make sure that the categories are mutually exclusive and cover everything, i.e. category total = gene total
    if not total_between_all_cases == len(gene_bin_counts):
        raise Exception("Something is wrong with bin counts! %s != %s (%s)"
                        %(total_between_all_cases, len(gene_bin_counts), len(bin_counts_simplified)))
    # TODO should probably unit-test all this!


def check_proper_enrichments(gene_stats_data, one_sided_thresholds=False, ratio_difference=3, min_alleles=2, 
                             FDR_cutoff=0.5, FDR_cutoff_2=0.1, ):
    """ Check to make sure the hits are enriched in the phenotype bins overall.

    Rather than e.g. depleted in the phenotype bin compared to wt, or enriched in one phenotype bin and depleted in another.

    For two-sided cases (where we're looking at both better-growing and worse-growing phenotypes), 
        check that one side has enrichment, and also check that BOTH sides don't have enrichment, because that would be weird.

    Only do this check for genes with FDR<cutoff; ignore cases where <min_alleles have a phenotype (on either side if two-sided).
    Consider something enriched if it's ratio_difference times larger than wt.  (Not super clear what this value should be...)
    """
    # ratios are slightly weird - we want 1/0 to be inf but 0/0 to be 0 (so 5:0:0 shows up as a hit on the 5:0 side but not the 0:0)
    def ratio_calc(phen_ratio,wt_ratio):
        if wt_ratio == 0:   return 0 if phen_ratio==0 else float('inf')
        else:               return phen_ratio/wt_ratio
    # Just doing two-sided and one-sided separately, since in the two-sided case we want to check both sides separately
    #   AND also make sure that we don't have enrichment on both sides.
    if one_sided_thresholds:
        bin_counts_simplified = {gene:(sum(d[0][:-1]), d[0][-1]) for (gene,d) in gene_stats_data.items()}
        bin_totals = [sum(x) for x in zip(*bin_counts_simplified.values())]
        for gene, bin_counts in bin_counts_simplified.items():
            FDR = gene_stats_data[gene][2]
            if FDR > FDR_cutoff or numpy.isnan(FDR): continue
            ratios = [ratio_calc(b,t) for b,t in zip(bin_counts,bin_totals)]
            ratio = ratio_calc(ratios[0], ratios[1])
            if ratio < ratio_difference:
                print ("Warning: gene %s (FDR %s, full bin counts %s, simplified %s:%s, "
                       +"phen:wt ratio normalized to totals %.2g:%.2g) doesn't look like a proper hit!")%(gene, FDR,
                                                              ':'.join(str(x) for x in gene_stats_data[gene][0]), 
                                                              bin_counts[0], bin_counts[1], ratios[0], ratios[1])
    else:
        N_bins = len(gene_stats_data.values()[0][0])
        middle_bin = int((N_bins-1)/2)
        bin_counts_simplified = {gene:(sum(d[0][:middle_bin]), d[0][middle_bin], sum(d[0][middle_bin+1:])) 
                                 for (gene,d) in gene_stats_data.items()}
        bin_totals = [sum(x) for x in zip(*bin_counts_simplified.values())]
        growth_count_data = [FDR_cutoff_2, 0, 0, 0, FDR_cutoff, 0, 0, 0]
        for gene, bin_counts in bin_counts_simplified.items():
            FDR = gene_stats_data[gene][2]
            if FDR > FDR_cutoff or numpy.isnan(FDR): continue
            ratios = [ratio_calc(b,t) for b,t in zip(bin_counts,bin_totals)]
            ratio1 = ratio_calc(ratios[0], ratios[1])
            ratio2 = ratio_calc(ratios[2], ratios[1])
            if max(ratio1, ratio2) < ratio_difference:
                print ("WARNING: gene %s (FDR %s, full bin counts %s, simplified %s, "
                       +"phen:wt ratios normalized to total %.2g/%.2g) doesn't look like a proper hit!")%(gene, FDR,
                                                              ':'.join(str(x) for x in gene_stats_data[gene][0]), 
                                                              ':'.join(str(x) for x in bin_counts), ratio1, ratio2)
            if min(ratio1, ratio2) > ratio_difference and min(bin_counts[0], bin_counts[2])>=min_alleles:
                print ("WARNING: gene %s (FDR %s, full bin counts %s, simplified %s, "
                       +"phen:wt ratios normalized to total %.2g/%.2g) "
                       +"looks like a hit with both better and worse growth!")%(gene, FDR,
                                                          ':'.join(str(x) for x in gene_stats_data[gene][0]), 
                                                          ':'.join(str(x) for x in bin_counts), ratio1, ratio2)
            # the better-vs-worse comparison can be done on straight numbers (ratios[0] etc) instead of ratios to wt (ratio1/2) -
            #  simpler that way! Otherwise in the ratio1==ratio2==inf case I'd have to do ratios[0]/ratios[2] comparisons anyway...
            if ratios[0] > ratios[2]:   growth_count_data[5] += 1
            elif ratios[0] < ratios[2]: growth_count_data[6] += 1
            else:                       growth_count_data[7] += 1
            if FDR < FDR_cutoff_2:
                if ratios[0] > ratios[2]:   growth_count_data[1] += 1
                elif ratios[0] < ratios[2]: growth_count_data[2] += 1
                else:                       growth_count_data[3] += 1
        print "Out of all putative hits: (FDR<%s) %s worse-growing, %s better-growing, %s mixed; (FDR<%s) %s, %s, %s."%tuple(
                                                                                                                growth_count_data)
    # TODO unit-test!


def scatterplot(readcountsS, readcountsC, sample_name, control_name, info='', same_scale=True, correlations=True):
    # I want min readcounts ABOVE zero so I can use them to set the symlog threshold
    min_readcount_x = min(x for x in readcountsC if x)
    min_readcount_y = min(x for x in readcountsS if x)
    if same_scale:
        min_readcount_x = min_readcount_y = min(min_readcount_x, min_readcount_y)
    mplt.figure(figsize=(8,8))
    mplt.scatter(readcountsC, readcountsS, edgecolor='None', color='k', s=20, alpha=0.2)
    mplt.xscale('symlog', linthreshx=min_readcount_x, linscalex=0.5)
    mplt.yscale('symlog', linthreshy=min_readcount_y, linscaley=0.5)
    scalemax_x = mplt.xlim()[1]
    scalemax_y = mplt.ylim()[1]
    if same_scale:
        scalemax_x = scalemax_y = max(scalemax_x, scalemax_y)
    mplt.xlim(-min_readcount_x/2, scalemax_x)
    mplt.ylim(-min_readcount_y/2, scalemax_y)
    # if one side IS all above zero (because of a readcount cutoff), use different x scale
    if not same_scale and min(readcountsC)>0:
        mplt.xlim(min_readcount_x/1.5, scalemax_x)
    mplt.xlabel('control %s readcount (log scale)'%control_name)
    mplt.ylabel('sample %s readcount (log scale)'%sample_name)
    if correlations:
        spearman_corr, spearman_pval = scipy.stats.spearmanr(readcountsS, readcountsC)
        pearson_corr, pearson_pval = scipy.stats.pearsonr(readcountsS, readcountsC)
        corr_info = '; correlation Spearman %.2f, Pearson %.2f'%(spearman_corr, pearson_corr)
    else:
        corr_info = ''
    mplt.title('%s vs %s readcount correlation, %s\n(each dot is an IB)%s'%(sample_name, control_name, info, corr_info))
    return min_readcount_x, min_readcount_y


def _plot_thresholds(phenotype_thresholds, linthreshx, linthreshy, color='dodgerblue', linestyle='dashed'):
    """ Plot phenotype thresholds on existing symlog plot.
    """
    xmax = mplt.xlim()[1]
    for ratio in phenotype_thresholds:
        if ratio==0 or numpy.isinf(ratio):  continue
        # if I plot a line on a symlog plot it'll be a straight line rather than switching between linear and log-scales correctly!
        # so to make these lines right I need to avoid putting anything in the linear region on either axis, 
        #   i.e. find the lowest x such that x >= linthreshx and ratio*x >= linthreshy.
        # IF I want the lines to just stop at the edge of the linear region and not do anything weird:
        #mplt.plot([edge_point, xmax], [edge_point*ratio, xmax*ratio], color=color, linestyle=linestyle)
        # BUT I also want to extend at least the outer lines down so that it's obvious what happens to x=0 and y=0 points
        #   even though they're out of the log-scale region...  To get the best approximation of that, 
        #   use 0, and the point where both the scales are linear, and the point where both are log.
        #   (this is still not perfect because in the region where one scale is linear and one is log the line should be a curve...)
        edge_points = [0, min(linthreshx, linthreshy/ratio), max(linthreshx, linthreshy/ratio), xmax]
        mplt.plot(edge_points, [ratio*x for x in edge_points], color=color, linestyle=linestyle)
        # But this is ugly... I could ALSO just fake it and make straight lines going all the way to -xmax, which looks better, 
        #   when there are no points in the affected regions anyway...  BUT sometimes there are points at y=0 that get crossed 
        #       by the lines even though they shouldn't be, so this isn't a good general solution!
        #   AND in the non-symmetrical case the lines are still not quite straight below edge_point! 
        #       Probably happens when the vertical and horizontal linear ranges are different sizes?  Not sure how to fix that... 
        #edge_point = max(linthreshx, linthreshy/ratio)
        #mplt.plot([-xmax, edge_point, xmax], [-xmax/ratio, edge_point*ratio, xmax*ratio], color=color, linestyle=linestyle)


def _plot_Poisson_thresholds(phenotype_thresholds, linthreshx, linthreshy, color='dodgerblue', linestyle='solid'):
    """ Plot phenotype thresholds on existing symlog plot.
    """
    for raw_offset in phenotype_thresholds:
        # for negative offsets just plot the absolute offset but on the other side
        offset = abs(raw_offset)
        # I don't want to plot every number up to 10k, so do 1,2,3,..., 100,200,300 etc, modified for the way I use offset.
        subrange = range(int(1*2**offset), int(100*2**offset), 1) + range(int(100*2**offset), int(mplt.xlim()[1]), 100)
        poisson_func = poisson_function(offset)
        poisson = [poisson_func(x) for x in subrange]
        if raw_offset>0:    mplt.plot(poisson, subrange, color=color, linestyle=linestyle)
        elif raw_offset<0:  mplt.plot(subrange, poisson, color=color, linestyle=linestyle)
        else:               print "PROBLEM: Weird offset value %s!"%raw_offset


def readcount_scatterplots(screen_data, screen_data_below_min, IBs_per_gene_filtered, gene_stats_data, outfolder, 
                           sample_name, control_name, phenotype_thresholds, min_reads, 
                           Poisson=False, symmetrical_min=False, min_reads_2=0, FDR_cutoffs=[0.5,0.05]):
    """ Make scatterplots of the sample-vs-control for all mutants (raw, norm, and hits), with cutoff lines drawn. 
    
    Three plots:
    1) raw readcounts with min readcount lines drawn, not filtered by gene/feature/confidence
    2) normalized readcounts filtered by min_readcount and gene/feature/etc, with phenotype threshold lines drawn
    3) same but only for multi-allele genes, highlighting all alleles of hit genes 
    Only plot filtered IBs (by library presence, feature, confidence, multiple-IBs-per-mutant). 
    """
    ### PLOT 1
    all_IBs_sorted = sorted(set(screen_data.keys() + screen_data_below_min.keys()))
    NaN = float('NaN')
    readcountsS = [screen_data.get(IB, screen_data_below_min.get(IB, NaN))[0] for IB in all_IBs_sorted]
    readcountsC = [screen_data.get(IB, screen_data_below_min.get(IB, NaN))[2] for IB in all_IBs_sorted]
    # MAYBE-TODO instead of showing entirely unfiltered IBs, maybe remove min_reads filtering but keep the rest?
    #   Honestly I should probably just have that be a separate function that happens before filtering in main...
    scatterplot(readcountsS, readcountsC, sample_name, control_name, 'raw')
    # if min_reads is a straight cutoff on one or both samples, plot lines, easy:
    if not symmetrical_min:
        mplt.plot((min_reads, min_reads), mplt.ylim(), color='dodgerblue', linestyle='dashed')
        if min_reads_2:
            mplt.plot(mplt.xlim(), (min_reads_2, min_reads_2), color='dodgerblue', linestyle='dashed')
    # if it's a sum of both, that's tricky to plot on a symlog plot, so just plot all points that add up to the min!
    #   and I have to use lists instead of just range generators because range generators fail to plot, ugh.
    else:
        mplt.plot([mplt.ylim()[0]] + list(range(min_reads+1)) + [min_reads], 
                  [min_reads] + list(reversed(range(min_reads+1))) + [mplt.xlim()[0]], color='dodgerblue', linestyle='dashed')
    plotting_utilities.savefig(outfolder + '/readcount-scatterplot__%s--%s__raw.png'%(sample_name, control_name))
    mplt.close()
    ### PLOT 2
    # TODO is there ANY sensible way to plot Poisson curves on normalized data??
    all_IBs_sorted = sorted(set.union(*IBs_per_gene_filtered.values()))
    readcountsS = [screen_data[IB][(0 if Poisson else 1)] for IB in all_IBs_sorted]
    readcountsC = [screen_data[IB][(2 if Poisson else 3)] for IB in all_IBs_sorted]
    linthreshx,linthreshy = scatterplot(readcountsS, readcountsC, sample_name, control_name, 
                                        ('raw' if Poisson else 'normalized'), same_scale=(symmetrical_min or Poisson))
    if Poisson: _plot_Poisson_thresholds(phenotype_thresholds, linthreshx, linthreshy)
    else:       _plot_thresholds(phenotype_thresholds, linthreshx, linthreshy)
    # if I want a diagonal (doesn't make sense on raw readcounts, but could help on normalized, if there are few thresholds):
    # _plot_thresholds([1], min(readcountsC), color='grey', linestyle='dotted')
    plotting_utilities.savefig(outfolder + '/readcount-scatterplot__%s--%s__thresholds.png'%(sample_name, control_name))
    # TODO TODO TODO implement optional Poisson color-coding by bin!
    ### PLOT 3
    if len(FDR_cutoffs) == 1:       colors, sizes = ['red'], [20]
    elif len(FDR_cutoffs) == 2:     colors, sizes = ['#aa0000', 'red'], [15, 30]  # green+limegreen is also decent. Colors are hard!
    else:                           raise Exception("I didn't define colors for %s FDR cutoffs for plot!"%len(FDR_cutoffs))
    for FDR_cutoff,color,size in zip(FDR_cutoffs,colors,sizes):
        hit_genes = set(gene for gene,data in gene_stats_data.items() if data[-1]<FDR_cutoff)
        if not hit_genes: continue
        hit_IBs_sorted = sorted(set.union(*[IBs for gene,IBs in IBs_per_gene_filtered.items() if gene in hit_genes]))
        readcountsS = [screen_data[IB][(0 if Poisson else 1)] for IB in hit_IBs_sorted]
        readcountsC = [screen_data[IB][(2 if Poisson else 3)] for IB in hit_IBs_sorted]
        mplt.scatter(readcountsC, readcountsS, edgecolor='None', color=color, s=size, alpha=1)
    mplt.title('%s vs %s readcount correlation with putative hit genes, '%(sample_name, control_name)
               +'\n(red = all alleles highlighted for genes under FDR cutoffs %s)'%(', '.join(str(x) for x in FDR_cutoffs))
               +'\n(grey/black = all alleles of genes with 2+ alleles only)')
    plotting_utilities.savefig(outfolder + '/readcount-scatterplot__%s--%s__hits.png'%(sample_name, control_name))
    mplt.close()
    # TODO plot a few alpha/size combinations maybe...
    # MAYBE-TODO  maybe also color-coding other high-phenotype alleles to show if they're in one-allele genes 
    #   or in genes with one hit allele and multiple non-hit alleles or what
    #  (need some phenotype threshold for alleles to count - can use same as in print_ratio_counts)


def gene_full_analysis(screen_sample_data, screen_control_data, library_data_by_IB, phenotype_thresholds, Poisson, min_reads, 
                       symmetrical_min=False, features="CDS,intron,5'UTR,5'UTR_intron", max_conf=4, min_alleles_for_FDR=1, 
                       sample_name='sample', control_name='control', scatterplot_outfolder=None, gene_names={}, 
                       one_sided_thresholds=False, min_reads_2=0, if_full_library=False, workspace_size=2e5):
    """ Do the basit statistical analysis as described in module docstring.

    Inputs:
     - screen_sample_data and screen_control_data - IB:(raw_readcount, normalized_readcount) dictionaries for the sample and control
     - library_data_by_IB - IB:insertion_namedtuple dict derived from the standard library rearray file
     - phenotype_thresholds - list of sample/control normalized ratio thresholds by which mutants will be binned
     - min_reads - the minimum number of reads in the control, below which mutants won't be included in the analysis to avoid noise
     - symmetrical_min - if True, min_reads applies to reads in control+sample, not just control
     - features - only mutants in one of those gene features will be included in the analysis
        (if a mutant is in multiple features of different splice variants, it's enough that one of the features is on the list)
     - max_conf - only mutants with a mapping confidence this or higher will be included in the analysis
     - min_alleles_for_FDR - p-values will be calculated for all genes, but FDRs will only be calculated for genes 
        with at least this many alleles, to avoid artificially increasing the FDRs by including genes that cannot be significant
     - min_reads_2 - the minimum number of reads in the SAMPLE, below which mutants won't be included in the analysis.
        This is normally not needed, and was just implemented to temporarily deal with a weird data situation.

    Output: gene:[binned_allele_counts, pval, FDR] dictionary.
    """
    all_IBs = set(IB for IB,x in screen_sample_data.items() if x[0]) | set(IB for IB,x in screen_control_data.items() if x[0])
    # screen_data is an IB:(raw_sample, norm_sample, raw_ctrl, norm_ctrl) readcount dict.
    screen_data = {IB: screen_sample_data.get(IB, (0,0)) + screen_control_data.get(IB, (0,0)) for IB in all_IBs}
    print_mutant_read_numbers(screen_data, sample_name, control_name)
    screen_data_below_min = filter_IBs_by_min_readcount(screen_data, min_reads, symmetrical_min, min_reads_2)
    if symmetrical_min: filter_text = 'Filtered (min %s control+sample reads)'%min_reads
    elif min_reads_2:   filter_text = 'Filtered (min %s control and %s sample reads)'%(min_reads, min_reads_2)
    else:               filter_text = 'Filtered (min %s control reads)'%min_reads
    print_mutant_read_numbers(screen_data, sample_name, control_name, filter_text)
    binned_IBs_by_phenotype = bin_IBs_by_phenotype(screen_data, phenotype_thresholds, Poisson)
    screen_lines_per_gene = raw_screen_data_per_gene(library_data_by_IB, screen_data, features, max_conf)
    IBs_per_gene_filtered = filter_IBs_per_gene(screen_lines_per_gene, symmetrical_min, library_data_by_IB, if_full_library)
    if not any(IBs_per_gene_filtered.values()):
        raise Exception("No IBs passed the filters!")
    gene_bin_counts = get_gene_bin_counts(IBs_per_gene_filtered, binned_IBs_by_phenotype)
    gene_stats_data = gene_statistics(gene_bin_counts, min_alleles_for_FDR, workspace_size)
    print "Alleles per gene (top 5, filtered): ", dict(collections.Counter(len(x) 
                                                                for x in IBs_per_gene_filtered.values()).most_common(5))
    print_ratio_counts(gene_bin_counts, .95, one_sided_thresholds)
    check_proper_enrichments(gene_stats_data, one_sided_thresholds, min_alleles=min_alleles_for_FDR)
    print "Number of hit genes by FDR cutoff:  " + ', '.join("%s: %s"%(x, sum(1 for d in gene_stats_data.values() if d[-1] <= x))
                                                              for x in (0.001, 0.01, 0.05, 0.1, 0.3, 0.5))
    if scatterplot_outfolder:
        readcount_scatterplots(screen_data, screen_data_below_min, IBs_per_gene_filtered, gene_stats_data, scatterplot_outfolder, 
                               sample_name, control_name, phenotype_thresholds, min_reads, Poisson, symmetrical_min, min_reads_2, 
                               FDR_cutoffs=[0.5,0.05])
    # sort the genes by FDR, then pval (FDR is NaN if <N alleles, so categorize those same as FDR=1, I guess)
    sorted_genes = sorted(gene_stats_data.items(), key = lambda (g,d): ((d[2] if not numpy.isnan(d[2]) else 1), d[1]))
    # print top 10 hits with decent FDRs, and if there are <10 of those, print one but not more with bad FDRs (bad=0.5).
    top10_hits = sorted_genes[:10]
    Nbad = sum(1 for g,d in top10_hits if not d[2]<=0.5)    # using not< instead of > to catch NaNs
    if Nbad>1:  top10_hits = top10_hits[:-Nbad+1]
    print "Top 10 hit genes (FDRs, bin counts): ", ', '.join(["%s (%.2g, %s)"%(gene_names.get(g, g), 
                                                                               d[2], ':'.join(str(x) for x in d[0])) 
                                                              for (g,d) in top10_hits])
    return gene_stats_data
    # TODO test!


def write_outfile(gene_stats_data, phenotype_thresholds, outfile, annotation, ann_header):
    """ write plaintext tab-separated output file containing all the genes and their binned allele numbers and p-values/FDRs.
    """
    header = ['gene']
    header += ['alleles_with_growth_%.2g-%.2g'%(a,b) for (a,b) in zip(phenotype_thresholds, phenotype_thresholds[1:])]
    header += 'raw_p-value FDR'.split()
    header += ann_header
    DATA = []
    for gene,(bin_counts,pval,FDR) in gene_stats_data.items():
        # make FDR a bit more readable
        if numpy.isnan(FDR):    FDR = 'NA'
        DATA.append([gene] + bin_counts + [pval, FDR] + annotation[gene])
    # sort by pval not FDR, so that the good-pval-no-FDR ones show up before pval=FDR=1 ones!  Otherwise the sort is the same anyway.
    DATA.sort(key = lambda x: (x[header.index('raw_p-value')], -x[1], x[0]))
    with open(outfile, 'w') as OUTFILE:
        OUTFILE.write('\t'.join(header) + '\n')
        for line in DATA:
            OUTFILE.write('\t'.join(str(x) for x in line) + '\n')
    return DATA, header


def compare_two_samples(sample1, control1, gene_stats_data1, sample2, control2, gene_stats_data2, outfile=None, minval=-1):
    """ Plot a scatterplot of gene FDRs between samples.

    Both sample/control should be names (strings). 
    Both gene_stats_data should be gene:(...,FDR) dicts (... is irrelevant), or pickle file names containing same.
    """
    if type(gene_stats_data1)==str:  gene_stats_data1 = general_utilities.unpickle(gene_stats_data1)
    if type(gene_stats_data2)==str:  gene_stats_data2 = general_utilities.unpickle(gene_stats_data2)
    all_genes_sorted = sorted(set(gene_stats_data1.keys()) | set(gene_stats_data2.keys()))
    FDRs1 = [-gene_stats_data1.get(gene,[2])[-1] for gene in all_genes_sorted]
    FDRs2 = [-gene_stats_data2.get(gene,[2])[-1] for gene in all_genes_sorted]
    FDRs1 = [-2 if numpy.isnan(x) else x for x in FDRs1]
    FDRs2 = [-2 if numpy.isnan(x) else x for x in FDRs2]
    mplt.figure(figsize=(8,8))
    mplt.scatter(FDRs1, FDRs2, edgecolor='None', color='k', s=50, alpha=0.3)
    mplt.xscale('symlog', linthreshx=0.001)
    mplt.yscale('symlog', linthreshy=0.001)
    mplt.xlim(minval-0.5, 0.0003)
    mplt.ylim(minval-0.5, 0.0003)
    mplt.xticks([0, -0.001, -0.01, -0.1, -1], '0 0.001 0.01 0.1 1'.split())
    mplt.yticks([0, -0.001, -0.01, -0.1, -1], '0 0.001 0.01 0.1 1'.split())
    mplt.title("False discovery rate (FDR) for %s vs %s (each dot is a gene)"%(sample1, sample2))
    mplt.xlabel('%s FDR (log scale)'%sample1)
    mplt.ylabel('%s FDR (log scale)'%sample2)
    for cutoff in (0.5, 0.05):
        mplt.plot(mplt.xlim(), (-cutoff, -cutoff), color='dodgerblue', linestyle='dashed')
        mplt.plot((-cutoff, -cutoff), mplt.ylim(), color='dodgerblue', linestyle='dashed')
    if outfile:
        plotting_utilities.savefig(outfile)
        mplt.close()


def main(args, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.

    This is written for very specific pickled input/output file formats mostly for my convenience:
        input file should be a pickled dictionary with sample names as keys and IB:readcount dictionaries as values.
        output will be the pickled output of gene_full_analysis.
    """
    # LATER-TODO if I want other people to use this, I should add plaintext input/output formats...
    try:
        [outfolder] = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: exactly one output folder is required!")
    try:                os.mkdir(outfolder)
    except OSError:     pass
    # LATER-TODO add options for library/annotation folder/filenames
    lib_folder = os.path.expanduser('~/experiments/arrayed_library/basic_library_data/')
    if options.full_library:
        library_table_header = general_utilities.unpickle(lib_folder+'large-lib_full_header.pickle')
        global Mutant5   # for some reason this has to be global or else the thing fails
        Mutant5 = collections.namedtuple('Mutant5', library_table_header)
        library_table = general_utilities.unpickle(lib_folder+'large-lib_full.pickle')
    else:
        library_table_header = general_utilities.unpickle(lib_folder+'large-lib_rearray_header.pickle')
        global Mutant   # for some reason this has to be global or else the thing fails
        Mutant = collections.namedtuple('Mutant', library_table_header)
        library_table = general_utilities.unpickle(lib_folder+'large-lib_rearray.pickle')
    library_data_by_IB = {x.IB: x for x in library_table}
    ann_file = os.path.expanduser('~/experiments/reference_data/chlamy_annotation/annotation+loc_data+header_v5.5.pickle')
    annotation, ann_header = general_utilities.unpickle(ann_file)
    gene_names = parse_annotation_file.best_gene_name_dict(annotation, ann_header)
    # run basic pipeline, with txt+pickle outfile; if replicates, analyze both and compare gene FDRs.
    phenotype_thresholds = modify_phenotype_thresholds(options.phenotype_thresholds, 
                                                       options.Poisson_thresholds, options.one_sided_thresholds)
    samples =  [x.strip() for x in options.sample.split(',')]
    controls = [x.strip() for x in options.control.split(',')]
    if len(samples) > 2:                        raise Exception("Can't have more than two samples!")
    if len(controls) > len(samples):            raise Exception("Can't have more controls than samples!")
    if len(samples)==2 and len(controls)==1:    controls *= 2
    if options.data_dictionary:
        screen_all_data = general_utilities.unpickle(options.data_dictionary)
    gene_stats_data_all = []
    for (sample, control) in zip(samples, controls):
        if options.data_dictionary:     screen_sample_data = screen_all_data[sample]
        else:                           screen_sample_data = general_utilities.unpickle(sample)
        if options.data_dictionary:     screen_control_data = screen_all_data[control]
        else:                           screen_control_data = general_utilities.unpickle(control)
        if sample.endswith('.pickle'):  sample = sample[:-len('.pickle')]
        if control.endswith('.pickle'): control = control[:-len('.pickle')]
        gene_stats_data = gene_full_analysis(screen_sample_data, screen_control_data, library_data_by_IB, 
                             phenotype_thresholds, options.Poisson_thresholds, 
                             options.min_reads, options.symmetrical_min_reads, options.features.split(','), options.max_conf, 
                             options.min_alleles_for_FDR, sample, control, outfolder, gene_names, 
                             options.one_sided_thresholds, options.min_reads_2, options.full_library, options.workspace_size)
        gene_stats_data_all.append((sample, control, gene_stats_data))
        # MAYBE-TODO change write_outfile to include multiple samples together?
        write_outfile(gene_stats_data, phenotype_thresholds, outfolder+'/results__%s--%s.txt'%(sample, control), 
                      annotation, ann_header)
    if len(samples) == 1:   general_utilities.pickle(gene_stats_data, outfolder+'/all_data.pickle')
    else:                   general_utilities.pickle(gene_stats_data_all, outfolder+'/all_data.pickle')
    # MAYBE-TODO maybe pickle the outfile/header too?  Together or separately from gene_stats_data
    # plot FDR comparison scatterplot for replicates (axes reversed so "good" cases are on top right
    if len(samples) == 2:
        compare_two_samples(*gene_stats_data_all[0]+gene_stats_data_all[1]+(outfolder + '/gene_FDR_scatterplot.png',))


def compare_multiples(gene_stats_datasets, names, overall_name='', outfile=None, linthresh=0.001, wraparound=True):
    """ Given a list of any number of gene results, plot the FDRs together: each gene is a line, each sample is a column)

    gene_stats_datasets should be either a list of gene:(...,FDR) dicts (... is irrelevant), or of pickle file names containing same.
    names should be a list of strings.
    """ 
    if len(set(type(x) for x in gene_stats_datasets)) != 1:
        raise Exception("Passed multiple types in the gene_stats_datasets list! %s"%(set(type(x) for x in gene_stats_datasets)))
    elif type(gene_stats_datasets[0])==str:
        gene_stats_datasets = [general_utilities.unpickle(x) for x in gene_stats_datasets]
    gene_all_FDRs = {}
    for gene in set.union(*[set(gene_stats_data.keys()) for gene_stats_data in gene_stats_datasets]):
        FDRs = [-gene_stats_data.get(gene,[2])[-1] for gene_stats_data in gene_stats_datasets]
        FDRs = [-2 if numpy.isnan(x) else x for x in FDRs]
        gene_all_FDRs[gene] = FDRs
    mplt.figure()
    for gene, FDRs in gene_all_FDRs.items():
        # plot -1 to N+1 and wrap the FDRs around, to see the extra connections and make all points have equal visual weight
        if not wraparound:
            mplt.plot(range(len(names)), FDRs, color='black', alpha=0.3)
        else:
            mplt.plot(range(-1, len(names)+1), [FDRs[-1]] + FDRs + [FDRs[0]], color='black', alpha=0.3)
    # LATER-TODO would be REALLY nice to color-code and label the top 10 genes...
    mplt.yscale('symlog', linthreshy=linthresh)     # the symlog step is very slow!
    mplt.ylim(-2, linthresh*0.3)
    mplt.xlim(-0.5, len(names)-0.5)
    mplt.yticks([0, -0.001, -0.01, -0.1, -1], '~0 0.001 0.01 0.1 1'.split())
    mplt.xticks(range(len(names)), names)
    mplt.title("FDR comparison for %s %s\n(each line is a gene, showing FDRs in all the samples)"%(overall_name, ', '.join(names)))
    mplt.ylabel('FDR (log scale; cutoff lines are 0.05 and 0.5)')
    mplt.xlabel('')
    for cutoff in (0.5, 0.05):
        mplt.plot(mplt.xlim(), (-cutoff, -cutoff), color='dodgerblue', linestyle='dashed')
    if outfile:
        plotting_utilities.savefig(outfile)
        mplt.close()


########################################### Additional analysis/display/plotting #################################################

def normalize_readcounts_single(sample_dict, norm_total=1e6):
    """ Given an IB:count dict, return matching one with normalized counts to a total of norm_total.
    """
    total = sum(sample_dict.values())
    return {IB: (count/total*norm_total) for (IB,count) in sample_dict.items()} 


def normalize_readcounts_multi(sample_dict, norm_total=1e6):
    """ Given a {sample:{IB:count}} dict, return matching one with normalized counts for each sample to a total of norm_total.
    """
    return {sample: normalize_readcounts_single(IB_counts, norm_total) for (sample, IB_counts) in sample_dict.items()}


### Notes and old interactive code for stuff I might want to implement later
def __NOTES():
    """
    ### more looking at IBs with a phenotype

    How many genes are they in?  Remember to split & cases and exclude intergenic/weirdness!
        >>> proper_gene = lambda IB: bool(FULL_by_IB[IB].feature != 'intergenic' 
        ...                               and FULL_by_IB[IB].gene not in mutant_IB_RISCC_classes.SPECIAL_GENE_CODES.all_codes)
        >>> [len(set.union(*[set(FULL_by_IB[IB].gene.split(' & ')) for IB in X if proper_gene(IB)])) 
        ...     for X in (SCREEN_IB_hits[(50, 0.1)])]
        [2821, 2565, 3245, 15350]

    Are any of those insertions confidence 5, or did we already exclude them?  Are any of them intergenic-only, cassette-only, etc?  What would the numbers be if we excluded those?
        >>> min_reads, ratio_cutoff = 50, 0.1
        >>> SCREEN_IBs_filtered = {IB for IB in SCREEN_data if FULL_by_IB[IB].confidence_level <= 4 and proper_gene(IB)}
        >>> hits1 = {IB for (IB,data) in SCREEN_data.items() if data[5][1]>=min_reads and data[3][1]/data[5][1]<=ratio_cutoff
        ...          and IB in SCREEN_IBs_filtered}
        >>> hits2 = {IB for (IB,data) in SCREEN_data.items() if data[5][1]>=min_reads and data[4][1]/data[5][1]<=ratio_cutoff
        ...          and IB in SCREEN_IBs_filtered}
        >>> hits12 = hits1 | hits2
        >>> all_IBs = {IB for (IB,data) in SCREEN_data.items() if data[5][1]>=min_reads}
        >>> print min_reads, ratio_cutoff, [(len(x), len(x)/len(SCREEN_IBs_filtered)) for x in (hits1, hits2, hits12, all_IBs)]
        50 0.1 [(3174, 0.030426780167951226), (2907, 0.027867249511100887), (3810, 0.036523639710111584), (192175, 1.8422389662180298)]

    How many mutants and genes are we looking at now?
        >>> [len(set((FULL_by_IB[IB].plate, FULL_by_IB[IB].well) for IB in X)) for X in (hits1, hits2, hits12, all_IBs)]
        [2638, 2369, 3109, 100874]
        >>> [len(set.union(*[set(FULL_by_IB[IB].gene.split(' & ')) for IB in X if proper_gene(IB)])) 
        ...     for X in (hits1, hits2, hits12, all_IBs)]
        [2249, 2046, 2599, 15350]

    # looking at overall hit/non-hit numbers - and why do we have so many hit mutants but so few hit genes?

    What are the hit/non-hit totals?  What's the % of hits?
        >>> gene_hits_nonhits_totals, gene_hits_total_fractions = {}, {}
        >>> for (cutoffs, data) in gene_hits_nonhits.items():
        ...     repl1 = [sum(h[0][x] for h in data.values()) for x in (0,1)]
        ...     repl2 = [sum(h[1][x] for h in data.values()) for x in (0,1)]
        ...     gene_hits_nonhits_totals[cutoffs] = (repl1, repl2)
        ...     gene_hits_total_fractions[cutoffs] = (repl1[0]/sum(repl1), repl2[0]/sum(repl2))
        ...     print cutoffs, repl1, repl1[0]/sum(repl1), repl2, repl2[0]/sum(repl2)
        (50, 0.1) [1656, 49425] 0.0324190990779 [1543, 49538] 0.0302069262544
        (30, 0.3) [5645, 45779] 0.109773646546 [5044, 46380] 0.0980864965775
        (20, 0.5) [11869, 39688] 0.23021122253 [10851, 40706] 0.210466086079

    Later I ran this again and the numbers were very slightly off (1655 instead of 1656 for 50, etc) - probably just some scipy change??

    For each version, look at a few things (just at replicate1, for simplicity):
     - how many genes have how many hits?
     - how many genes have X hits and 0 non-hits?
     - how many genes have X hits and Y non-hits with more hits than non?
        >>> for (cutoffs, data) in gene_hits_nonhits.items():
        ...     print cutoffs
        ...     print collections.Counter(x[0][0] for x in data.values())
        ...     print collections.Counter(x[0][0] for x in data.values() if x[0][1]==0)
        ...     print collections.Counter((x[0][0],x[0][1]) for x in data.values() if x[0][1]<x[0][0])
        * 50, 0.1
        {0: 11135, 1: 1181, 2: 120, 3: 29, 4: 8, 5: 6, 7: 2, 9: 2, 10: 1, 15: 1, 29: 1}
        {1: 129, 0: 83, 2: 14, 3: 6, 5: 1}
        {(1, 0): 129, (2, 0): 14, (2, 1): 13, (3, 0): 6, (3, 1): 2, (5, 1): 2, (10, 8): 1, (3, 2): 1, (4, 3): 1, (5, 0): 1, (4, 2): 1, (4, 1): 1}
        * 30, 0.3
        {0: 8635, 1: 2872, 2: 681, 3: 172, 4: 53, 5: 34, 6: 11, 7: 6, 9: 6, 8: 4, 11: 2, 10: 1, 12: 1, 13: 1, 15: 1, 108: 1, 17: 1, 52: 1, 21: 1, 20: 1, 29: 1}
        {1: 393, 0: 57, 2: 37, 3: 9, 5: 1, 6: 1}
        {(1, 0): 393, (2, 1): 59, (2, 0): 37, (3, 1): 12, (3, 2): 10, (3, 0): 9, (4, 2): 5, (4, 3): 3, (4, 1): 3, (9, 3): 2, (8, 6): 1, (6, 0): 1, (17, 1): 1, (5, 0): 1, (5, 1): 1, (5, 3): 1}
        * 20, 0.5
        {0: 6053, 1: 3936, 2: 1418, 3: 579, 4: 229, 5: 130, 6: 50, 7: 28, 8: 10, 9: 8, 14: 8, 11: 6, 12: 5, 10: 4, 13: 2, 15: 2, 18: 2, 20: 2, 16: 1, 17: 1, 21: 1, 22: 1, 25: 1, 26: 1, 28: 1, 29: 1, 39: 1, 42: 1, 47: 1, 54: 1, 248: 1, 122: 1}
        {1: 827, 2: 126, 0: 43, 3: 33, 5: 3, 4: 2, 6: 1, 18: 1}
        {(1, 0): 827, (2, 1): 207, (2, 0): 126, (3, 2): 60, (3, 1): 38, (3, 0): 33, (4, 3): 22, (4, 2): 15, (4, 1): 13, (5, 4): 11, (5, 3): 9, (6, 5): 6, (5, 2): 6, (5, 1): 5, (5, 0): 3, (6, 4): 3, (4, 0): 2, (7, 6): 2, (8, 6): 2, (18, 0): 1, (6, 2): 1, (6, 3): 1, (13, 11): 1, (9, 3): 1, (6, 0): 1, (7, 4): 1, (10, 2): 1}

    How many genes have ANY hit IBs in either replicate?
        >>> sum(1 for x in gene_hits_nonhits[(50, 0.1)].values() if x[0][0] or x[1][0])
        1601
    But earlier I calculated that all the IBs with phenotypes cover 2599 genes - why is this number different?  Probably because we excluded 3'UTRs here and included them earlier.
    Trying the previous method again but excluding 3'UTRs (roughly, because to do it properly we'd have to split the gene+feature for "X & Y" cases - so we're excluding slightly too many?):
        >>> len(set.union(*[set(library_data_by_IB[IB].gene.split(' & ')) for IB in hits12 if proper_gene(IB) and "3" not in library_data_by_IB[IB].feature])) 
        1868
    Okay, that's still more than 1601 - LATER-TODO why?  Confidence level 5 is already filtered out here...

    Anyway, even disregarding that problem, how come we have thousands of mutants with phenotypes but only end up with 44 hit genes?  There are 1601 genes with at least one insertion with a phenotype, 3'UTRs excluded.  How many insertions in each replicate?
        >>> for repl in (0,1):
        ...     genes = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0])
        ...     alleles = sum(x[repl][0] for x in gene_hits_nonhits[(50,0.1)].values())
        ...     genes_with_2 = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0] >= 2)
        ...     with_1 = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0] == 1)
        ...     with_morehits = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0] >= x[repl][1])
        ...     with_2_and_morehits = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0] >= max(x[repl][1],2))
        ...     with_3 = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0] >= 3)
        ...     with_1_0 = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0] == 1 and x[repl][1] == 0)
        ...     with_1_1 = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0] == 1 and x[repl][1] == 1)
        ...     with_1_more = sum(1 for x in gene_hits_nonhits[(50,0.1)].values() if x[repl][0] == 1 and x[repl][1] > 1)
        ...     print genes, alleles, genes_with_2, with_1, with_morehits, with_2_and_morehits, with_3, with_1_0, with_1_1, with_1_more
        1350 1655 170 1180 416 54 50 130 149 901
        1264 1542 155 1109 394 51 45 126 134 849

    So basically >1k genes have an allele with a phenotype, but only <200 of genes have TWO alleles with phenotypes, and most of those have more non-phenotype than phenotype alleles; only ~50 genes have THREE alleles with phenotypes.  So we have a lot of cases of either genes with a single allele, or genes that have many alleles but only one of them is a hit, likely due to second-site mutations OR wrong mapping.  It looks like ~900 genes are in the second category (1 allele with a phenotype, 2+ without), ~140 are 1-vs-1, ~130 1-vs-0.  

    So the main problem is all the genes with only a single allele with a phenotype and many alleles without a phenotype.  Why might that happen?
    1) phenotype of the one allele could be due to second-site mutation
    2) phenotype of the one allele could be due to that allele being wrongly mapped to that gene
    3) the gene could be a real hit with a mild phenotype, so most of its alleles would be below our threshold
    4) anything else?
    MAYBE-TODO try to figure out which case it is?  (Probably some of all three, but at least which one is dominant)
    1) what's the confidence level distribution of the insertions that are the single alleles with a phenotype for a gene?
    2) do more of them have a mapped second insertion, compared to the overall distribution of #insertions per mutant?
        If they do have a mapped second insertion, how often is it in a hit or candidate gene, or known photosynthesis gene, etc?
    3) do the other insertions in those genes have a mild phenotype on average, or are they wildtype-like?

    Plot hits/nonhits scatterplots, with statistically significant genes on top in green:
        >>> for (cutoffs, data_both_replicates) in GENE_stats.items():
        ...  for (N,data) in enumerate(data_both_replicates):
        ...     info = (N+1,) + cutoffs
        ...     mplt.figure(figsize=(8,8))
        ...     hits, nonhits, pvals = [[data[g][x] for g in genes_sorted] for x in (0,1,-1)]
        ...     mplt.scatter(hits, nonhits, marker='.', c='k', edgecolors='None', s=300, alpha=0.2)
        ...     mplt.scatter([h for (h,p) in zip(hits,pvals) if p<0.05], 
        ...                  [h for (h,p) in zip(nonhits,pvals) if p<0.05], 
        ...                  marker='.', c='limegreen', edgecolors='None', s=300, alpha=1)
        ...     mplt.xscale('symlog')
        ...     mplt.yscale('symlog')
        ...     scalemax = max(max(hits), max(nonhits))*2
        ...     mplt.xlim(-1, scalemax)
        ...     mplt.ylim(-1, scalemax)
        ...     mplt.xlabel('number of hit alleles per gene')
        ...     mplt.ylabel('number of non-hit alleles per gene')
        ...     mplt.title('Gene hit/non-hit alleles (green=statistically significant hit gene)' 
        ...                +'\nreplicate %s, min %s reads, ratio cutoff %s'%info)
        ...     plotting_utilities.savefig('Xiaobo_screen/gene-hits-nonhits_repl%s_%s-%s.png'%info)
        ...     mplt.close()

    Plot correlation between the two replicate p-values:
        >>> for (cutoffs, data) in GENE_stats.items():
        ...     mplt.figure(figsize=(8,8))
        ...     mplt.scatter([data[0][g][-1] for g in genes_sorted], 
        ...                  [data[1][g][-1] for g in genes_sorted], 
        ...                  marker='.', c='k', edgecolors='None', s=300, alpha=0.2)
        ...     mplt.xscale('log')
        ...     mplt.yscale('log')
        ...     mplt.xlim(1e-7, 2)
        ...     mplt.ylim(1e-7, 2)
        ...     mplt.xlabel('replicate 1 gene p-values')
        ...     mplt.ylabel('replicate 2 gene p-values')
        ...     mplt.title('Comparing FDR-adjusted p-values between replicates, min %s reads, ratio cutoff %s'%cutoffs)
        ...     plotting_utilities.savefig('Xiaobo_screen/gene-pval_repl-correlation_%s-%s.png'%cutoffs)
        ...     mplt.close()

    What are the outliers, and do we have any idea why they look like that?
        >>> [(g,x) for (g,x) in GENE_stats[(50,0.1)][0].items() if x[3]<0.03 and x[3]/GENE_stats[(50,0.1)][1][g][3] < 0.3]
        [('Cre13.g581850', (5, 5, 7.9363976597733572e-06, 0.014156265882847163)), ('Cre02.g111550', (10, 8, 4.5549322566182101e-11, 5.687288415613497e-07))]

    LATER-TODO do any mutants in any of the hits have insertions in another hit?  If so, the phenotype might be due to that and one of the hits may not be real.

    MAYBE-TODO try redoing all this with using flagellar proteome instead of all-genes as the comparison set?  Maybe only take mutants that have a flagellar proteome insertion and NO other insertions - OR not, because maybe that biases the comparison against mutants with multiple insertions, which are more likely to show a phenotype?


    ################### looking at the final gene statistics, picking out/printing hits

    What are the totals for 1 allele and 2+ alleles? (This uses the old GENE_stats which gets overwritten next)
        >>> print len(genes_sorted_1all), len(genes_sorted_2all), len(gene_hits_nonhits_pvals_1), len(GENE_stats[c][0])
        3542 8851 12393 12476
    Okay, so most genes had 2+ alleles, but ~30% had 1.  But why is the last number bigger??  The old analysis had more genes total??
        >>> print len(data), len(genes_sorted), sum(1 for g in genes_sorted if sum(data[g][0]) >= 1)
        12476 12476 12393
    HA, so we were even including genes with zero alleles??  Apparently yes, a few!  Oops.
        >>> print sum(1 for g in genes_sorted if sum(data[g][0]) >= 1 or sum(data[g][1]) >= 1)
        12393
        >>> print sum(1 for g in genes_sorted if sum(data[g][0]) == 0 and sum(data[g][1]) == 0)
        83

    Grab the top hits, with 0.3 pval cutoff this time - should be same set of genes
        >>> HITS_old = HITS_set
        >>> GENE_stats_final = (gene_hits_nonhits_pvals_1, gene_hits_nonhits_pvals_2)
        >>> hit_pval_cutoff = 0.3
        >>> HITS1, HITS2 = [{g:x for (g,x) in GENE_stats_final[i].items() if x[-1]<=hit_pval_cutoff} for i in (0,1)]
        >>> HITS_set = set(HITS1.keys()) | set(HITS2.keys())
        >>> HITS_12_data = {g: (GENE_stats_final[0][g], GENE_stats_final[1][g]) for g in HITS_set}
        >>> print len(HITS1), len(HITS2), len(HITS_12_data), len(set(HITS1.keys()) & set(HITS2.keys())), HITS_set == HITS_old
        37 34 44 27 True
    Yes, same set of genes, good!

    How many alleles do these mostly have?  Median is 3, often 2, sometimes a lot more.
        >>> collections.Counter([x[0][0]+x[0][1] for x in HITS_12_data.values()])
        Counter({2: 18, 3: 6, 6: 5, 5: 4, 4: 3, 7: 3, 10: 2, 12: 1, 14: 1, 18: 1})
        >>> numpy.median([x[0][0]+x[0][1] for x in HITS_12_data.values()])
        3.0

    How many total mutants with a phenotype do they have?  Or total mutants overall?
        >>> print sum([x[0][0] for x in HITS_12_data.values()]), sum([x[0][0]+x[0][1] for x in HITS_12_data.values()])
        127 201

    Print the hits nicely:
        >>> for (g, ((h1, n1, p1, q1), (h2, n2, p2, q2))) in sorted(HITS_12_data.items(), key = lambda (x,(y,z)): min(y[-1],z[-1])):
        ...     n, d, D = annotation_data[g][:3]
        ...     a = annotation_data[g][-1]
        ...     print "%s   replicate1 %s:%s pval %.3g   replicate2 %s:%s pval %.3g (%s; %s; %s; %s)"%(g, h1, n1, q1, h2, n2, q2, 
        ...         n, d, D, a)
    Cre07.g316050   replicate1 2:0 pval 0.364      replicate2 1:1 pval 1 (CDJ2; Chloroplast DnaJ-like protein)
    Cre10.g430150   replicate1 2:0 pval 0.364      replicate2 1:1 pval 1 (-; -; -; tetratricopeptide repeat-containing protein)
    Cre10.g429650   replicate1 2:0 pval 0.364      replicate2 1:1 pval 1 (-; -; -; -)

    And save them to a file for Xiaobo:
        >>> with open('Xiaobo_screen/new_hits_renormalized.txt', 'w') as OUTFILE:
        ...     OUTFILE.write('gene repl1_hits repl1_non-hits repl1_pval repl2_hits repl2_non-hits repl2_pval'.replace(' ','\t')
        ...                   + '\t' + '\t'.join(annotation_header) + '\n')
        ...     for (g, ((h1, n1, p1, q1), (h2, n2, p2, q2))) in sorted(HITS_12_data.items(), key = lambda (x,(y,z)): min(y[-1],z[-1])):
        ...         A = annotation_data[g]
        ...         OUTFILE.write('\t'.join([str(x) for x in (g, h1, n1, q1, h2, n2, q2)] + A) + '\n')

    Also print file with all ALLELES of all these genes, and including their readcounts as well as their rearray information!  Rearray only, for now, for people who want to pick them (NOTE that not all will really still be growing!).  Sort by GENE and then allele.
        >>> all_hit_IBs = set(x.IB for x in REARRAY_filtered if strip_v55(x.gene) in HITS_set)
        >>> all_hit_alleles_full = []
        >>> for IB in all_hit_IBs:
        ...     if IB not in SCREEN_data:   continue
        ...     A = REARRAY_by_IB[IB]
        ...     line = [A.mutant_ID, A.strand, A.side, A.chromosome, A.min_position, A.full_position, A.gene, A.orientation, A.feature, 
        ...              A.confidence_level, A.LEAPseq_distance, A.LEAPseq_percent_confirming, A.IB]
        ...     for B in SCREEN_data[IB][:-1]:
        ...         line += [B[0], B[1]]
        ...     all_hit_alleles_full.append(line)
        ... 
        >>> all_hit_alleles_full.sort(key = lambda x: (x[6], x[1], x[3], x[4]))
        >>> with open('Xiaobo_screen/all_hit_alleles_rearray.txt', 'w') as OUTFILE:
        ...     OUTFILE.write('\t'.join(cpl3_header[:10] + ['LEAPseq_distance', 'LEAPseq_percent_confirming'] + cpl3_header[10:]) + '\n')
        ...     for line in all_hit_alleles_full:
        ...         OUTFILE.write('\t'.join(str(x) for x in line) + '\n')

    MAYBE-TODO redo everything below with new FDR-adjusted p-values if needed!  I redid some things but not everything; some things don't really matter.

    # generating more liberal candidate/potential hit list for further testing by crossing

    For a larger list for crossing, Xiaobo picked all genes with a p-value lower than a [1,1] one (different for each set), so replicate that too (just for the best version).  Look at raw p-val and FDR-adjusted q-val:
        >>> [set((p,q) for (g, (h,n,p,q)) in GENE_stats_final[i].items() if h==n==1) for i in (0,1)]
        [set([(0.063896078782841745, 1.0)]), set([(0.059603075366071533, 1.0)])]
    That q-value is 1.0 in both cases, so I'll have to use the raw p-value as a cutoff:
    AND also check that the ratio of hits:nonhits is actually higher than expected, not lower! Otherwise some results will be for lower ratios, which is silly.
        >>> maybeHITS1 = {g:x for (g,x) in GENE_stats_final[0].items() if x[-2]<0.063896078782841745
        ...              and (x[1]==0 or x[0]/x[1] > gene_hits_nonhits_totals[(50,0.1)][0][0]/gene_hits_nonhits_totals[(50,0.1)][0][1])}
        >>> maybeHITS2 = {g:x for (g,x) in GENE_stats_final[1].items() if x[-2]<0.059603075366071533
        ...              and (x[1]==0 or x[0]/x[1] > gene_hits_nonhits_totals[(50,0.1)][1][0]/gene_hits_nonhits_totals[(50,0.1)][1][1])}
        >>> maybeHITS_set = set(maybeHITS1.keys()) | set(maybeHITS2.keys())
        >>> maybeHITS_12_data = {g:(GENE_stats_final[0][g], GENE_stats_final[1][g]) for g in maybeHITS_set}
        >>> print len(maybeHITS1), len(maybeHITS2), len(maybeHITS_12_data)
        247 230 308

    What are the numbers without the "real" hits?
        >>> print len(maybeHITS1) - len(HITS1), len(maybeHITS2) - len(HITS2), len(maybeHITS_set - HITS_set)
        210 196 264

    What are the highest p-values here?  Okay, good, so we can use 0.058 as a general cutoff despite differences between replicates.
        >>> max(x1[-2] for (x1,x2) in maybeHITS_12_data.values() if x1[-2] < 0.063896078782841745)
        0.056147563292997332
        >>> max(x2[-2] for (x1,x2) in maybeHITS_12_data.values() if x2[-2] < 0.059603075366071533)
        0.057267084988197589

    How many just 1-vs-0 in either replicate?
        >>> only_1_0_r1 = set(g for (g,x) in GENE_stats_final[0].items() if x[0]==1 and x[1]==0)
        >>> only_1_0_r2 = set(g for (g,x) in GENE_stats_final[1].items() if x[0]==1 and x[1]==0)
        >>> print len(only_1_0_r1), len(only_1_0_r2), len(only_1_0_r1 | only_1_0_r2), len(only_1_0_r1 & only_1_0_r2)

    Pickle this:
        >>> general_utilities.pickle((maybeHITS_12_data, maybeHITS_set, maybeHITS1, maybeHITS2), 'pickled_data/maybeHITS.pickle')
    Unpickle if needed:
        >>> (maybeHITS_12_data, maybeHITS_set, maybeHITS1, maybeHITS2) = general_utilities.unpickle('pickled_data/maybeHITS.pickle')

    Print a list of those for Xiaobo - same format as main hits:
        >>> with open('Xiaobo_screen/possible_hits_renormalized.txt', 'w') as OUTFILE:
        ...     OUTFILE.write('gene repl1_hits repl1_non-hits repl1_pval repl2_hits repl2_non-hits repl2_pval'.replace(' ','\t')
        ...                   + '\t' + '\t'.join(annotation_header) + '\n')
        ...     for (g, ((h1, n1, p1, q1), (h2, n2, p2, q2))) in sorted(maybeHITS_12_data.items(), 
        ...                                 key = lambda (x,(y,z)): min(y[-1],z[-1])):
        ...         A = annotation_data[g]
        ...         OUTFILE.write('\t'.join([str(x) for x in (g, h1, n1, q1, h2, n2, q2)] + A) + '\n')

    """
    pass

########################################### Tests #################################################


class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__(self):
        sys.exit("NO UNIT-TESTS FOR THIS MODULE")
    # LATER-TODO add unit-tests!


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    from testing_utilities import run_functional_tests
    test_folder = "test_data"
    sys.exit("NO TESTS DEFINED!")
    # tests in (testname, [test_description,] arg_and_infile_string) format
    test_runs = [ ]
    # argument_converter converts (parser,options,args) to the correct argument order for main
    argument_converter = lambda parser,options,args: (args, options)
    # use my custom function to run all the tests, auto-detect reference files, compare to output.
    return run_functional_tests(test_runs, define_option_parser(), main, test_folder, 
                                argument_converter=argument_converter, append_to_outfilenames='.txt') 
    # LATER-TODO add run-tests!


if __name__=='__main__':
    parser = define_option_parser()
    options,args = parser.parse_args()

    # if run with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        #print("\n * unit-tests for the ______ module")
        #test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(______
        #unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        print("\n * unit-tests for this module (%s)"%sys.argv[0])
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise pass the arguments to the main function
    main(args, options)
