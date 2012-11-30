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
import basic_seq_utilities

DEFAULT_NUCLEAR_GENOME_FILE = '~/experiments/reference_data/genomes_and_indexes/Chlre4-nm.fa'
DEFAULT_ALL_GENOME_FILE = '~/experiments/reference_data/genomes_and_indexes/Chlre4-nm_chl-mit.fa'
DEFAULT_GENOME_CASSETTE_FILE = '~/experiments/reference_data/genomes_and_indexes/Chlre4-nm_chl-mit_cassette-pMJ013b.fa'

########################################### Mutant simulation functions ##########################################

### mappability over genome

# TODO!
# TODO should we look at 21bp chunks, or 20bp, or a combination of both?  Do we have any mutants with only 20bp flanking regions? Check!


### simulate dataset with N randomly positioned mutants, taking into account mappable/unmappable positions!

# TODO


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

   
def gene_counts_for_mutant_subsets(dataset, step_size=100, max_N_mutants=3):
    """ Return numbers of genes with N mutants for different-sized random subsets of dataset.

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list of mutants (mutant_analysis_classes.Insertional_mutant instances).

    Return a N:list_of_gene_numbers dictionary, where N is each value between 1 and max_N_mutants, 
     and list_of_gene_numbers contains the number of genes with at least N mutants 
      in randomly chosen subsets of dataset, starting at 0 and going up to the full dataset size in step_size steps.
    (So if dataset contains 200 mutants, and step_size is 100, each list will have 3 values, for 0, 100 and 200 mutants.)
    (Note that the last value is not for the full dataset size, but for the closest lower number divisible by step_size.)
    """
    # extract list of mutants from dataset; randomize the order
    if isinstance(dataset, mutant_analysis_classes.Insertional_mutant_pool_dataset):    mutants = list(dataset)
    else:                                                                               mutants = dataset
    random.shuffle(mutants)
    # get the gene counts for each mutant number
    gene_counts = {}
    for N_mutants in range(1,max_N_mutants+1):
        gene_counts[N_mutants] = [_genes_with_N_mutants(mutants[:subset_size], N_mutants) 
                                  for subset_size in range(0, len(mutants), step_size)]
    return gene_counts


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
