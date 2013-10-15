#! /usr/bin/env python2.7

"""
VERY TEMPORARY PROGRAM to make joint tab-separated growthrate files (with annotation).
 -- Weronika Patena, 2012
"""

# TODO really I should just rewrite growth rates to be a sensible part of mutant data using subclassing, and then mutant_join_datasets.py can do this.

# standard library
from __future__ import division
import sys
import unittest
# other packages
from numpy import mean, median, isnan, isinf, isneginf
# my modules
import mutant_growth_rates
import mutant_analysis_classes

def all_nan_growthrates(mutants):
    for mutant in mutants:
        try:  mutant.growth_rate
        except AttributeError:  mutant.growth_rate = float('nan')
    
if __name__=='__main__':

    if len(sys.argv) < 6:
        sys.exit("ERROR: needs at least five arguments!  Two or more input files, a mutant file containing all the datasets pooled, an annotation file, and an output file name.")

    all_infiles = sys.argv[1:-3]
    reffile = sys.argv[-3]
    annfile = sys.argv[-2]
    outfile = sys.argv[-1]

    print "Infiles %s, outfile %s"%(all_infiles, outfile)

    # read all infile datasets, set missing growthrates to NaN
    all_datasets = [mutant_growth_rates._read_growthrate_data_from_file_tmp(infile) for infile in all_infiles]
    for dataset in all_datasets:
        all_nan_growthrates(dataset)
    
    # grab all positions so we can iterate over them later
    all_mutant_positions = sorted(reduce( set.union, [set([m.position for m in dataset]) for dataset in all_datasets], set([]) ))
    print "total %s mutants found"%len(all_mutant_positions)

    # Now grab the readcount file with all the mutants it in, so we can get the gene info etc from that!  And add gene-info annotation to it.
    mutant_data = mutant_analysis_classes.Insertional_mutant_pool_dataset(infile=reffile)
    mutant_data.add_gene_annotation(annfile, True)

    # print data to outfile
    with open(outfile, 'w') as OUTFILE:
        # print header: mutant info (no full_position), then growthrate info (RZ-34 min100, min300, min1000, RZ-67 same), annotation.  I made the headers by hand here - TODO that's awful, fix!
        mutant_info_header = "# chromosome\tstrand\tmin_position\tgene\torientation\tfeature\tmain_sequence"  
        growthrate_header = "growthrate_RZ-34_min100\tgrowthrate_RZ-34_min300\tgrowthrate_RZ-34_min1000\tgrowthrate_RZ-67_min100\tgrowthrate_RZ-67_min300\tgrowthrate_RZ-67_min1000"
        annotation_header = "PFAM\tPanther\tKOG\tKEGG_ec\tKEGG_Orthology\tbest_arabidopsis_TAIR10_hit_name\tbest_arabidopsis_TAIR10_hit_symbol\tbest_arabidopsis_TAIR10_hit_defline"
        OUTFILE.write("%s\t%s\t%s\n"%(mutant_info_header, growthrate_header, annotation_header))

        # For each position that was found in the dataset, print all the info to OUTFILE!  
        # Make all growth rates 0 if they're lower
        for position in all_mutant_positions:
            mutant = mutant_data.get_mutant(position)
            OUTFILE.write('\t'.join([position.chromosome, position.strand, str(position.min_position), mutant.gene, mutant.orientation, mutant.gene_feature, mutant.get_main_sequence()[0]]))
            OUTFILE.write('\t')
            for dataset in all_datasets:
                try:  
                    g = dataset.get_mutant(position).growth_rate
                except AttributeError:  
                    g = float('nan')
                if g < 0:
                    g = 0
                OUTFILE.write('%s\t'%g)
            # add gene annotation if present, or the made-up empty annotation if not
            try:
                if mutant.gene_annotation:
                    OUTFILE.write('\t'.join(mutant.gene_annotation))
                else:
                    OUTFILE.write('NO GENE DATA\t\t\t\t\t\t\t')
            except AttributeError:
                OUTFILE.write('NO GENE DATA\t\t\t\t\t\t\t')
            OUTFILE.write('\n')

