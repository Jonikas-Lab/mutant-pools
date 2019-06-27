#! /usr/bin/env python2.7

"""
Basic analysis pipeline for pooled mutant screens, resulting in p-values for each gene based on allele numbers and phenotypes: 
    - filter IBs by minimum control readcount; calculate sample:control phenotype ratios for each IB, 
    - sort mutants into bins based on the ratios and cutoffs (non-hits, weak/medium/strong hits)
    - for each gene, do a Fisher's exact test comparing the numbers of mutants in this gene in each bin to all mutants, 
        optionally filtering by gene feature and confidence level
    - do FDR-correction on the resulting p-values, optionally excluding genes with <N alleles

Input file should be a python pickle file containing a {sample_name:{IB:(raw_readcount, normalized_readcount)}} dictionary.
Output file is also a python pickle file.

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
import numpy
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
    parser.add_option('-p', '--phenotype_thresholds', type='string', default='0.0625,0.125,0.25,0.5', metavar='X,Y', 
              help="List of sample/control ratio thresholds for bins, comma-separated, no spaces, lowest first (default %default)."
                  +"The real thresholds used will be made two-sided, i.e. '10' -> '0.1,10' etc, unless -P is used.")
    parser.add_option('-P', '--one_sided_thresholds', action='store_true', default=False,
              help="Only look for slower-growing and not faster-growing mutants (default: two-sided)")
    parser.add_option('-m', '--min_reads', type='int', default='50', metavar='M', 
              help="Minimum raw control reads to include IB in analysis (default %default).")
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
    return parser


########################################### Basic pipeline #################################################


def filter_IBs_by_min_readcount(screen_data, min_readcount, min_readcount_2=0):
    """ Modifies screen_data to remove any IBs below the min readcounts.
    """
    for (IB,data) in screen_data.items():
        if data[2] < min_readcount or data[0] < min_readcount_2:
            del screen_data[IB]


def bin_IBs_by_phenotype(screen_data, phenotype_thresholds):
    binned_IBs = [set(IB for IB,data in screen_data.items() if min(a,b) <= data[1]/data[3] < max(a,b)) 
                  for (a,b) in zip(phenotype_thresholds, phenotype_thresholds[1:])]
    print "numbers of IBs in each phenotype bin: ", [len(x) for x in binned_IBs]
    return binned_IBs


def raw_screen_data_per_gene(library_data_by_IB, screen_data, features, max_conf):
    """ Make a list of screen per-IB data lines for each gene, filtering by feature and confidence
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


def filter_IBs_per_gene(screen_lines_per_gene, library_data_by_IB, if_full_library=False):
    """ Filter lists of screen data for each gene to remove multiple ones in the same mutant (keep the highest-control-readcount one)
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
            line_set.sort(key = lambda x: x[1][2])
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


def gene_statistics(gene_bin_counts, min_alleles_for_FDR=1):
    """ Calculate a p-value (using Fisher's exact test) and FDR (BH method) for each gene compared to all alleles

    Output: gene:[binned_allele_counts, pval, FDR] dictionary, where FDR is NaN for genes with <min_alleles_for_FDR alleles.
    """
    bin_totals = [sum(x) for x in zip(*gene_bin_counts.values())]
    print "numbers of filtered IBs in each phenotype bin: ", bin_totals
    # if there are only 2 bins, can use scipy.stats.fisher_exact; otherwise use my custom one that goes through R
    if len(bin_totals) == 1:    fisher_exact = lambda x: scipy.stats.fisher_exact(x)[1]
    else:                       fisher_exact = statistics_utilities.fisher_exact
    gene_pvals = {g:fisher_exact([bin_counts,bin_totals]) for (g,bin_counts) in gene_bin_counts.items()}
    if any(numpy.isnan(x) for x in gene_pvals.values()): raise Exception("Some pvals are NaN! Look into this!")
    # the FDR-correction has to be done on a list of pvalues, so separate out the genes that meet min_alleles_for_FDR
    genes_with_enough_alleles = sorted(g for (g,bin_counts) in gene_bin_counts.items() if sum(bin_counts) > min_alleles_for_FDR)
    FDRs_for_some_genes = statistics_utilities.FDR_adjust_pvalues([gene_pvals[g] for g in genes_with_enough_alleles], method='BH')
    FDR_dict = collections.defaultdict(lambda: float('NaN'), zip(genes_with_enough_alleles, FDRs_for_some_genes))
    # join the bin counts, pvals and FDRs into a single dictionary, with NaN FDRs as needed
    gene_stats_data = {g: [bin_counts, gene_pvals[g], FDR_dict[g]] for (g,bin_counts) in gene_bin_counts.items()}
    return gene_stats_data


def print_mutant_read_numbers(screen_data, header='Samples'):
    sumF = lambda x: general_utilities.format_number(sum(x), 0, 1)
    lenF = lambda x: general_utilities.format_number(len(x), 0, 1)
    overlap_IBs = set(IB for IB,x in screen_data.items() if x[0] and x[2])
    print "%s: %s %s mutants (%s reads);  Control %s %s (%s);  Overlap %s (%s/%s)"%(header,
       options.sample_key,  sumF(1 for x in screen_data.values() if x[0]), sumF(x[0] for x in screen_data.values()),
       options.control_key, sumF(1 for x in screen_data.values() if x[2]), sumF(x[2] for x in screen_data.values()),
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
    for (name, min_wt, max_wt, min_ph, max_ph) in [('3+:0',   3, inf,   0, 0),
                                                   ('2:0',    2, 2,     0, 0),   
                                                   ('1:0',    1, 1,     0, 0),   
                                                   ('3+:1',   3, inf,   1, 1),   
                                                   ('2:1',    2, 2,     1, 1),   
                                                   ('1:1',    1, 1,     1, 1),   
                                                   ('2+:2+',  2, inf,   2, inf), 
                                                   ('1:2+',   1, 1,     2, inf), 
                                                   ('0:1',    0, 0,     1, 1),   
                                                   ('0:2+',   0, 0,     2, inf)]: 
        N = sum(1 for wt,ph in bin_counts_simplified if min_wt<=wt<=max_wt and min_ph<=ph<=max_ph)
        print_data.append("%s: %s"%(name, N))
        total_between_all_cases += N
    print "Genes by wt:phenotype allele counts: " + ', '.join(print_data)
    # make sure that the categories are mutually exclusive and cover everything, i.e. category total = gene total
    if not total_between_all_cases == len(gene_bin_counts):
        raise Exception("Something is wrong with bin counts! %s != %s (%s)"
                        %(total_between_all_cases, len(gene_bin_counts), len(bin_counts_simplified)))
    # TODO should probably unit-test all this!


def check_proper_enrichments(gene_stats_data, one_sided_thresholds=False, min_alleles=2, FDR_cutoff=0.5, ratio_difference=3):
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
            ratios = [b/t for b,t in zip(bin_counts,bin_totals)]
            ratio = ratio_calc(ratios[0], ratios[1])
            if ratio < ratio_difference:
                print ratios, ratio, ratio_difference, ratio<ratio_difference
                print ("Warning: gene %s (FDR %s, full bin counts %s, simplified %s:%s, "
                       +"phenotype:wt ratio normalized to totals %.2g:%.2g) doesn't look like a proper hit!")%(gene, FDR,
                                                              ':'.join(str(x) for x in gene_stats_data[gene][0]), 
                                                              bin_counts[0], bin_counts[1], ratios[0], ratios[1])
    else:
        N_bins = len(gene_stats_data.values()[0][0])
        middle_bin = int((N_bins-1)/2)
        bin_counts_simplified = {gene:(sum(d[0][:middle_bin]), d[0][middle_bin], sum(d[0][middle_bin+1:])) 
                                 for (gene,d) in gene_stats_data.items()}
        bin_totals = [sum(x) for x in zip(*bin_counts_simplified.values())]
        worse_growing, better_growing = 0, 0
        for gene, bin_counts in bin_counts_simplified.items():
            FDR = gene_stats_data[gene][2]
            if FDR > FDR_cutoff or numpy.isnan(FDR): continue
            ratios = [b/t for b,t in zip(bin_counts,bin_totals)]
            ratio1 = ratio_calc(ratios[0], ratios[1])
            ratio2 = ratio_calc(ratios[2], ratios[1])
            if max(ratio1, ratio2) < ratio_difference:
                print ("WARNING: gene %s (FDR %s, full bin counts %s, simplified %s, "
                       +"phenotype:wt ratios normalized to total %.2g/%.2g) doesn't look like a proper hit!")%(gene, FDR,
                                                              ':'.join(str(x) for x in gene_stats_data[gene][0]), 
                                                              ':'.join(str(x) for x in bin_counts), ratio1, ratio2)
            if min(ratio1, ratio2) > ratio_difference and min(bin_counts[0], bin_counts[2])>=min_alleles:
                print ("WARNING: gene %s (FDR %s, full bin counts %s, simplified %s, "
                       +"phenotype:wt ratios normalized to total %.2g/%.2g) "
                       +"looks like a hit with both better and worse growth!")%(gene, FDR,
                                                          ':'.join(str(x) for x in gene_stats_data[gene][0]), 
                                                          ':'.join(str(x) for x in bin_counts), ratio1, ratio2)
            if ratio1 > ratio2:     worse_growing += 1
            elif ratio1 < ratio2:   better_growing += 1
        print "Out of all putative hits: %s worse-growing, %s better-growing (FDR<%s)"%(worse_growing, better_growing, FDR_cutoff)
    # TODO unit-test!


def gene_full_analysis(screen_sample_data, screen_control_data, library_data_by_IB, phenotype_thresholds, min_reads, features, 
               max_conf, gene_names={}, min_alleles_for_FDR=1, one_sided_thresholds=False, min_reads_2=0, if_full_library=False):
    """ Do the basit statistical analysis as described in module docstring.

    Inputs:
     - screen_sample_data and screen_control_data - IB:(raw_readcount, normalized_readcount) dictionaries for the sample and control
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
    if min(phenotype_thresholds) <= 1 <= max(phenotype_thresholds):
        raise Exception("All thresholds should be >1 or all <1! (current values are %s). "%phenotype_thresholds
                        +"Symmetrical ones will be added automatically to make the test two-sided if one_sided_thresholds is False.")
    # in one-sided cases, always make it so the wt bin is at the end - makes display better and everything else easier.
    if one_sided_thresholds:
        if numpy.mean(phenotype_thresholds) < 1:
            phenotype_thresholds = [0] + sorted(phenotype_thresholds) + [float('inf')]
        else:
            phenotype_thresholds = [float('inf')] + sorted(phenotype_thresholds, reverse=True) + [0]
    else:
        phenotype_thresholds = sorted([0] + [1/x for x in phenotype_thresholds] + phenotype_thresholds + [float('inf')])
    print "phenotype thresholds: ", phenotype_thresholds
    all_IBs =     set(IB for IB,x in screen_sample_data.items() if x[0]) | set(IB for IB,x in screen_control_data.items() if x[0])
    screen_data = {IB: screen_sample_data.get(IB, (0,0)) + screen_control_data.get(IB, (0,0)) for IB in all_IBs}
    print_mutant_read_numbers(screen_data)
    filter_IBs_by_min_readcount(screen_data, min_reads, min_reads_2)
    print_mutant_read_numbers(screen_data, 'Filtered (min %s control reads)'%min_reads)
    binned_IBs_by_phenotype = bin_IBs_by_phenotype(screen_data, phenotype_thresholds)
    screen_lines_per_gene = raw_screen_data_per_gene(library_data_by_IB, screen_data, features, max_conf)
    IBs_per_gene_filtered = filter_IBs_per_gene(screen_lines_per_gene, library_data_by_IB, if_full_library)
    gene_bin_counts = get_gene_bin_counts(IBs_per_gene_filtered, binned_IBs_by_phenotype)
    gene_stats_data = gene_statistics(gene_bin_counts, min_alleles_for_FDR)
    print "Alleles per gene (top 5, filtered): ", dict(collections.Counter(len(x) 
                                                                for x in IBs_per_gene_filtered.values()).most_common(5))
    print_ratio_counts(gene_bin_counts, .95, one_sided_thresholds)
    check_proper_enrichments(gene_stats_data, one_sided_thresholds, min_alleles=min_alleles_for_FDR)
    print "Number of hit genes by FDR cutoff:  " + ', '.join("%s: %s"%(x, sum(1 for d in gene_stats_data.values() if d[-1] <= x))
                                                              for x in (0.001, 0.01, 0.05, 0.1, 0.3, 0.5))
    # sort the genes by FDR, then pval (FDR is NaN if <N alleles, so categorize those same as FDR=1, I guess)
    sorted_genes = sorted(gene_stats_data.items(), key = lambda (g,d): ((d[2] if not numpy.isnan(d[2]) else 1), d[1]))
    print "Top 5 hit genes (FDRs, bin counts): ", ', '.join(["%s (%.2g, %s)"%(gene_names.get(g, g), 
                                                                              d[2], ':'.join(str(x) for x in d[0])) 
                                                             for (g,d) in sorted_genes[:5] if d[2]<1])
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
    ann_file = os.path.expanduser('~/experiments/reference_data/chlamy_annotation/annotation_data_and_header_v5.5.pickle')
    annotation, ann_header = general_utilities.unpickle(ann_file)
    gene_names = {gene:(a[0] if a[0]!='-' else (a[11] if a[11]!='-' else gene)) for gene,a in annotation.items()}
    gene_names = {g:n.split(',')[0] for g,n in gene_names.items()}
    # run basic pipeline, pickle output to outfile
    phenotype_thresholds = [float(x) for x in options.phenotype_thresholds.split(',')]
    phenotype_thresholds = [int(x) if round(x)==x else x for x in phenotype_thresholds]
    gene_stats_data = gene_full_analysis(screen_sample_data, screen_control_data, library_data_by_IB, phenotype_thresholds, 
                                         options.min_reads, options.features.split(','), 
                                         options.max_conf, gene_names, options.min_alleles_for_FDR, 
                                         options.one_sided_thresholds, options.min_reads_2, options.full_library)
    # TODO write output file
    # TODO make scatterplot with the readcount/cutoff lines drawn?
    general_utilities.pickle(gene_stats_data, outfile)


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
