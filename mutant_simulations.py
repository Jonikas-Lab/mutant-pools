#! /usr/bin/env python
"""
Plotting utilities specifically for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import unittest
import math
import random
from collections import defaultdict
# other packages
import numpy
import scipy
# my modules
import mutant_analysis_classes
import general_utilities
import seq_basic_utilities


########################################### Mutant simulation functions ##########################################

### mappability over genome

# TODO!


### simulate dataset with N randomly positioned mutants, taking into account mappable/unmappable positions!

# TODO


### simulate dataset with N randomly positioned mutants, matching the gap-size distribution of the real dataset, and taking into account mappable/unmappable positions

# LATER-TODO
# MAYBE-TODO would it be possible to also match the real hotspots and coldspots of the real dataset?


### number of genes (with one, two, etc mutants) vs number of mutants (randomly chosen mutant subsets)

def gene_counts_for_mutant_subsets(dataset, step_size, repeat_N):
    """ ____ """
    # TODO implement!
    pass


### TODO all this starting here until the TESTING section is just code copied from another file! Rewrite it to make sense and be modular; some of it should probably go in mutant_plotting_utilities.py

total_mutants = len(dataset)
total_genes = dataset.summary.total_genes_in_genome

def genes_with_N_mutants(mutants, N):
    gene_mutant_counts = defaultdict(int)
    for m in mutants:
        if m.gene not in mutant_analysis_classes.SPECIAL_GENE_CODES.all_codes:
            gene_mutant_counts[m.gene] += 1
    return len([1 for count in gene_mutant_counts.values() if count>=N])
   
def genes_with_no_mutants(mutants):
   return total_genes - genes_with_N_mutants(mutants, 1)
   

def make_random_data():
  random.shuffle(dataset._mutants_by_position.values())
  gene_counts = {}
  gene_counts[1] = [genes_with_N_mutants(mutants[:x], 1) for x in range(0, len(mutants), 100)]
  gene_counts[2] = [genes_with_N_mutants(mutants[:x], 2) for x in range(0, len(mutants), 100)]
  gene_counts[3] = [genes_with_N_mutants(mutants[:x], 3) for x in range(0, len(mutants), 100)]
  return gene_counts


colors = {1: 'm', 2: 'c', 3: 'y'}

# Plot it once with labels and make a legend:
gene_counts_by_Nmutants = make_random_data()
for N_mutants,gene_counts in gene_counts_by_Nmutants.items():
  mplt.plot(gene_counts, '.', linewidth=0, c=colors[N_mutants], label = "genes with %s mutants"%N_mutants)

mplt.legend(loc=2)

# Plot it all again 10 times with new random mutant subsets, to make sure we have a good coverage of the random space 
for _ in range(10):
  gene_counts_by_Nmutants = make_random_data()
  for N_mutants,gene_counts in gene_counts_by_Nmutants.items():
    mplt.plot(gene_counts, '.', linewidth=0, c=colors[N_mutants])


mplt.title('Number of genes hit vs number of mutants sequenced\n(chose random mutant subsets ten times, plotted all - they overlap)')
mplt.ylabel("Number of genes hit (out of %s total chlamy nuclear genes)"%total_genes)
mplt.yticks(mplt.yticks()[0], ["%i (%.0f%%)"%(x, x*100/total_genes) for x in mplt.yticks()[0]])
mplt.xlabel("Number of mutants (randomly chosen out of %s total)\n(counting only unique-genomic mutants - 36%% of total reads, 49%% of \"good\" reads)"%total_mutants)
mplt.xticks(mplt.xticks()[0], ["%i"%(x*100) for x in mplt.xticks()[0]])
mplt.ylim(-100, 7500)
mplt.xlim(-2, 140)
savefig('genes-hit_vs_mutant-number.png')


################################################# Testing etc ##################################################

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__(self):
        sys.exit("NO UNIT-TESTS FOR THIS MODULE")
    # LATER-TODO add unit-tests!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
