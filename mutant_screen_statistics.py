#! /usr/bin/env python2.7

"""
Basic analysis pipeline for pooled mutant screens, resulting in p-values for each gene based on allele numbers and phenotypes: 
    - normalize IB readcounts, calculate ratios between samples for each IB, 
    - sort mutants into bins based on the ratios (non-hits, weak/medium/strong hits) with filtering based on a raw readcount minimum
    - for each gene, do a Fisher's exact test comparing the numbers of mutants in this gene in each bin to all mutants, 
        optionally filtering by gene feature and confidence level
    - do FDR-correction on the resulting p-values, optionally excluding genes with <N alleles
Currently this module contains functions to do all these things, but each function has to be separately run in the interactive python shell (see ../../1802_Moshe+Frieder_new-screens/notes.txt for examples - the functions were first used there).

The code is largely based on a simpler analysis of earlier screen data done for Xiaobo in ../../../arrayed_library/1603_large-lib_figures-1/notes.txt section "### do statistics to find potential gene hits"; the more complex analysis method (with multiple thresholds for Fisher's exact test) is based on the analysis Robert did for the large-scale screen paper (2018?) - see code in ../../1802_Moshe+Frieder_new-screens/Robert_pipeline, although his code isn't very readable for me and I didn't use it.

 -- Weronika Patena, 2018

USAGE: mutant_screen_statistics.py [options] infile outfile
"""

# standard library
from __future__ import division
import sys
import os
import unittest
import collections
# other packages
import scipy.stats
# my modules
import general_utilities
import mutant_utilities
import mutant_IB_RISCC_classes
import statistics_utilities


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
    parser.add_option('-S', '--sample_key', type='string', default='', metavar='"Sample 1"', 
      help="Name of the sample in the infile dictionary (required) - put quotes around it if it has spaces or weird characters.")
    parser.add_option('-C', '--control_key', type='string', default='', metavar='"Sample 2"', 
      help="Name of the control in the infile dictionary (required) - put quotes around it if it has spaces or weird characters.")
    parser.add_option('-p', '--phenotype_thresholds', type='string', default='0.1,10', metavar='X,Y', 
              help="List of sample/control ratio thresholds for bins, comma-separated, no spaces, lowest first (default %default).")
    parser.add_option('-m', '--min_reads', type='int', default='50', metavar='M', 
              help="Minimum raw control reads to include IB in analysis (default %default).")
    parser.add_option('-x', '--min_reads_2', type='int', default='0', metavar='M', 
              help="Minimum raw SAMPLE reads to include IB in analysis - rarely needed! (default %default).")
    parser.add_option('-f', '--features', type='string', default="CDS,intron,5'UTR,5'UTR_intron", metavar='F,F',
              help="Only mutants in these gene features will be included in the analysis (default %default)")
    parser.add_option('-c', '--max_conf', type='int', default=4, metavar='C',
              help="Only mutants with mapping confidence from 1 to C will be included in the analysis (default %default)")
    parser.add_option('-a', '--min_alleles_for_FDR', type='int', default=1, metavar='A',
              help="Only genes with at least A qualified alleles are included in FDR correction (default %default)")
    return parser


########################################### Basic pipeline #################################################


def normalize_readcounts_single(sample_dict, norm_total=1e6):
    """ Given an IB:count dict, return matching one with normalized counts to a total of norm_total.
    """
    total = sum(sample_dict.values())
    return {IB: (count/total*norm_total) for (IB,count) in sample_dict.items()} 


def normalize_readcounts_multi(sample_dict, norm_total=1e6):
    """ Given a {sample:{IB:count}} dict, return matching one with normalized counts for each sample to a total of norm_total.
    """
    return {sample: normalize_readcounts_single(IB_counts, norm_total) for (sample, IB_counts) in sample_dict.items()}


def _bin_IBs_by_phenotype(screen_data, phenotype_thresholds, min_readcount, min_readcount_2=0):
    phenotype_thresholds = [0] + sorted(phenotype_thresholds) + [float('inf')]
    print "phenotype thresholds: ", phenotype_thresholds
    binned_IBs = [set(IB for IB,data in screen_data.items() if data[2] >= min_readcount and data[0] >= min_readcount_2
                      and low <= data[1]/data[3] < high) 
                  for (low, high) in zip(phenotype_thresholds, phenotype_thresholds[1:])]
    print "numbers of IBs in each phenotype bin: ", [len(x) for x in binned_IBs]
    return binned_IBs


def _raw_screen_data_per_gene(library_data_by_IB, screen_data, features, max_conf):
    """ Make a list of screen per-IB data lines for each gene, filtering by feature and confidence
    """
    screen_lines_per_gene = collections.defaultdict(list)
    for IB in screen_data.keys():
        x = library_data_by_IB[IB]
        if x.gene not in mutant_IB_RISCC_classes.SPECIAL_GENE_CODES.all_codes and x.feature != 'intergenic':
          if x.confidence_level <= max_conf:
            for (gene,feature) in zip(x.gene.split(' & '), x.feature.split(' & ')):
              if features is None or mutant_utilities.if_right_feature(feature, features):
                gene = mutant_utilities.strip_version(gene)
                screen_lines_per_gene[mutant_utilities.strip_version(gene)].append((IB,screen_data[IB]))
    return screen_lines_per_gene


def _filter_IBs_per_gene(screen_lines_per_gene, library_data_by_IB):
    """ Filter lists of screen data for each gene to remove multiple ones in the same mutant (keep the highest-control-readcount one)
    """
    IBs_per_gene_filtered = {}
    for gene,screen_lines in screen_lines_per_gene.items():
        screen_lines_by_mutant = collections.defaultdict(list)
        for (IB,data) in screen_lines:
            screen_lines_by_mutant[library_data_by_IB[IB].mutant_ID].append((IB,data))
        filtered_IBs = set()
        for line_set in screen_lines_by_mutant.values():
            line_set.sort(key = lambda x: x[1][2])
            filtered_IBs.add(line_set[-1][0])
        IBs_per_gene_filtered[gene] = filtered_IBs
    print "top 10 numbers of filtered alleles per gene: ", dict(collections.Counter(len(x) 
                                                                for x in IBs_per_gene_filtered.values()).most_common(10))
    return IBs_per_gene_filtered


def _get_gene_bin_counts(IBs_per_gene_filtered, binned_IBs_by_phenotype):
    """ Get counts of alleles in each hit/nonhit bin for each gene

    Genes with no alleles in any bins aren't included.
    """
    gene_bin_counts = {}
    for (gene, filtered_IBs) in IBs_per_gene_filtered.items():
        bin_counts = [len(filtered_IBs & bin_IBs) for bin_IBs in binned_IBs_by_phenotype]
        if sum(bin_counts):
            gene_bin_counts[gene] = bin_counts
    return gene_bin_counts


def _gene_statistics(gene_bin_counts, min_alleles_for_FDR=1):
    """ Calculate a p-value (using Fisher's exact test) and FDR (BH method) for each gene compared to all alleles
    """
    bin_totals = [sum(x) for x in zip(*gene_bin_counts.values())]
    print "numbers of filtered IBs in each phenotype bin: ", bin_totals
    # if there are only 2 bins, can use scipy.stats.fisher_exact; otherwise use my custom one that goes through R
    if len(bin_totals) == 1:    fisher_exact = lambda x: scipy.stats.fisher_exact(x)[1]
    else:                       fisher_exact = statistics_utilities.fisher_exact
    gene_pvals = {g:fisher_exact([bin_counts,bin_totals]) for (g,bin_counts) in gene_bin_counts.items()}
    # the FDR-correction has to be done on a list of pvalues, so separate out the genes that meet min_alleles_for_FDR
    genes_with_enough_alleles = sorted(g for (g,bin_counts) in gene_bin_counts.items() if sum(bin_counts) > min_alleles_for_FDR)
    FDRs_for_some_genes = statistics_utilities.FDR_adjust_pvalues([gene_pvals[g] for g in genes_with_enough_alleles], method='BH')
    FDR_dict = collections.defaultdict(lambda: float('NaN'), zip(genes_with_enough_alleles, FDRs_for_some_genes))
    # join the bin counts, pvals and FDRs into a single dictionary, with NaN FDRs as needed
    gene_stats_data = {g: [bin_counts, gene_pvals[g], FDR_dict[g]] for (g,bin_counts) in gene_bin_counts.items()}
    return gene_stats_data


def gene_full_analysis(screen_sample_data, screen_control_data, library_data_by_IB, 
                       phenotype_thresholds, min_reads, features, max_conf, min_alleles_for_FDR=1, min_reads_2=0):
    """ Do the basit statistical analysis as described in module docstring.

    Inputs:
     - screen_sample_data and screen_control_data - IB:raw_readcount dictionaries for the sample and control
     - library_data_by_IB - IB:insertion_namedtuple dict derived from the standard library rearray file
     - phenotype_thresholds - list of sample/control normalized ratio thresholds by which mutants will be binned
     - min_reads - the minimum number of reads in the control, below which mutants won't be included in the analysis to avoid noise
     - features - only mutants in one of those gene features will be included in the analysis
        (if a mutant is in multiple features of different splice variants, it's enough that one of the features is on the list)
     - max_conf - only mutants with a mapping confidence this or higher will be included in the analysis
     - min_alleles_for_FDR - p-values will be calculated for all genes, but FDRs will only be calculated for genes 
        with at least this many alleles, to avoid artificially increasing the FDRs by including genes that cannot be significant
     - min_reads_2 - the minimum number of reads in the SAMPLE, below which mutants won't be included in the analysis.
        This is normally not needed, and was just implemented to temporarily deal with a weird data situation.

    Output: gene:[binned_allele_counts, pval, FDR] dictionary.
    """
    screen_sample_data_norm = normalize_readcounts_single(screen_sample_data)
    screen_control_data_norm = normalize_readcounts_single(screen_control_data)
    all_IBs = set(screen_sample_data.keys()) | set(screen_control_data.keys())
    screen_data = {IB: [x.get(IB, 0) for x in (screen_sample_data, screen_sample_data_norm, 
                                               screen_control_data, screen_control_data_norm)] for IB in all_IBs}
    binned_IBs_by_phenotype = _bin_IBs_by_phenotype(screen_data, phenotype_thresholds, min_reads, min_reads_2)
    screen_lines_per_gene = _raw_screen_data_per_gene(library_data_by_IB, screen_data, features, max_conf)
    IBs_per_gene_filtered = _filter_IBs_per_gene(screen_lines_per_gene, library_data_by_IB)
    gene_bin_counts = _get_gene_bin_counts(IBs_per_gene_filtered, binned_IBs_by_phenotype)
    gene_stats_data = _gene_statistics(gene_bin_counts, min_alleles_for_FDR)
    print "number of hit genes by FDR cutoff: ", {x: sum(1 for d in gene_stats_data.values() if d[-1] <= x) 
                                                  for x in (0.3, 0.1, 0.05, 0.01, 0.001)}
    print "top 5 hit genes (with FDRs): ", ' '.join(["%s (%.2g)"%(g,d[-1]) for (g,d) in sorted(gene_stats_data.items(), 
                                                   key = lambda (g,d): d[-1] if not scipy.isnan(d[-1]) else 2) if d[-1] < 1] [:5])
    # MAYBE-TODO add annotation to get gene names?
    return gene_stats_data
    # TODO test!


def main(args, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.

    This is written for very specific pickled input/output file formats mostly for my convenience:
        input file should be a pickled dictionary with sample names as keys and IB:readcount dictionaries as values.
        output will be the pickled output of gene_full_analysis.
    """
    # LATER-TODO if I want other people to use this, I should add plaintext input/output formats...
    try:
        infile, outfile = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: exactly one infile and outfile are required!")
    screen_all_data = general_utilities.unpickle(infile)
    screen_sample_data = screen_all_data[options.sample_key]
    screen_control_data = screen_all_data[options.control_key]
    overlap_IBs = set(IB for (IB,x) in screen_sample_data.items() if x) & set(IB for (IB,x) in screen_control_data.items() if x)
    print "Sample %s: %s mutants, %s reads;  Control %s: %s mutants, %s reads"%(
            options.sample_key,  sum(x for x in screen_sample_data.values() if x),  sum(screen_sample_data.values()), 
            options.control_key, sum(x for x in screen_control_data.values() if x), sum(screen_control_data.values())) 
    print "  Overlap %s mutants, %s/%s reads"%(len(overlap_IBs), sum(screen_sample_data[x]  for x in overlap_IBs), 
                                                                 sum(screen_control_data[x] for x in overlap_IBs))
    # LATER-TODO add options for library folder/filenames
    lib_folder = os.path.expanduser('~/experiments/arrayed_library/basic_library_data/')
    REARRAY_table_header = general_utilities.unpickle(lib_folder+'large-lib_rearray_header.pickle')
    global Mutant   # for some reason this has to be global or else the thing fails
    Mutant = collections.namedtuple('Mutant', REARRAY_table_header)
    REARRAY_table = general_utilities.unpickle(lib_folder+'large-lib_rearray.pickle')
    library_data_by_IB = {x.IB: x for x in REARRAY_table}
    # run basic pipeline, pickle output to outfile
    gene_stats_data = gene_full_analysis(screen_sample_data, screen_control_data, library_data_by_IB, 
                                         phenotype_thresholds=[float(x) for x in options.phenotype_thresholds.split(',')], 
                                         min_reads=options.min_reads, features=options.features.split(','), 
                                         max_conf=options.max_conf, min_alleles_for_FDR=options.min_alleles_for_FDR,
                                         min_reads_2=options.min_reads_2)
    general_utilities.pickle(gene_stats_data, outfile)


########################################### Additional analysis/display/plotting #################################################

### Notes and old interactive code for stuff I might want to implement later
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