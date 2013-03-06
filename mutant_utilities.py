#! /usr/bin/env python

"""
Various help functions for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import sys
import unittest
import os
from collections import defaultdict
# other packages
import numpy
# my modules
import general_utilities
import basic_seq_utilities
import mutant_analysis_classes

DEFAULT_NUCLEAR_GENOME_FILE = os.path.expanduser('~/experiments/reference_data/genomes_and_indexes/Chlre4_notmasked.fa')
DEFAULT_ALL_GENOME_FILE = os.path.expanduser('~/experiments/reference_data/genomes_and_indexes/Chlre4-nm_chl-mit.fa')
DEFAULT_GENOME_CASSETTE_FILE = os.path.expanduser('~/experiments/reference_data/genomes_and_indexes/Chlre4-nm_chl-mit_cassette-pMJ013b.fa')
DEFAULT_GENE_POS_FILE = os.path.expanduser('~/experiments/reference_data/chlamy_annotation/Creinhardtii_169_gene.gff3')
DEFAULT_GENE_ANNOTATION_FILE = os.path.expanduser('~/experiments/reference_data/chlamy_annotation/Creinhardtii_169_annotation_info.txt')

STRAND_VAR_VALUES = ('+', '-', 'both', None)
DEFAULT_BIN_SIZE = 20000


def get_chromosome_lengths(genome_file=None):
    """ Return chromosome:length dictionary based on reading a genome fasta file. """
    original_input = genome_file
    if genome_file is None:
        genome_file = DEFAULT_GENOME_CASSETTE_FILE
    chromosome_lengths = defaultdict(int)
    try:
        for header,seq in basic_seq_utilities.parse_fasta(genome_file):
            chromosome_lengths[header] = len(seq)
        return dict(chromosome_lengths)
    except IOError:
        file_info = "default " if original_input is None else ""
        raise ValueError("%sgenome fasta file $s not found! Provide filename."%(file_info, genome_file))
    # MAYBE-TODO should this be in basic_seq_utilities or somewhere?  Except for the specific default value...


def get_mutant_positions_from_dataset(dataset, strand=None):
    """  Return chromosome_name:mutant_position_list for dataset.

    Dataset must be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list/set/something of mutant_analysis_classes.Insertional_mutant instances.
    Use the (known or assumed) position of the base before the insertion (min_position).

    If strand is None, take mutants regardless of strand; 
     if it's '+', '-' or 'both', take only mutants on that strand (both-stranded mutants are opposite-tandems); 
     all other strand values are illegal.
    """
    if not strand in STRAND_VAR_VALUES:
        raise ValueError("Illegal strand value %s! Must be one of %s"%(strand, STRAND_VAR_VALUES))
    chromosome_position_dict = defaultdict(list)
    for mutant in dataset:
        position = mutant.position
        if strand is None or position.strand==strand:
            chromosome_position_dict[position.chromosome].append(position.min_position)
    return chromosome_position_dict


def get_histogram_data_from_positions(position_dict, bin_size=DEFAULT_BIN_SIZE, chromosome_lengths=None, chromosomes=None, 
                                      first_bin_offset=0, special_last_bin=True, merge_last_bin_cutoff=0.5, normalize_last_bin=True):
    """ Given a chromosome:position_list dict, return a chromosome:counts_per_bin dict, with counts_per_bin a numpy array. 
    
    The input can actually be either a chromosome:position_list or a (chromosome,strand):position_list dict - 
     the results will be the same for either, the different-strand lists will just be added together. 

    The positions in position_list will be binned into bin_size-sized bins over the length of each chromosome, 
     giving counts_per_bin numpy array, with length matching the number of bins in the chromosome. 
    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default file will be used) - use the lengths to decide how many bins there will be in the chromosome.

    Skip and ignore the first first_bin_offset of each chromosome.
    If special_last_bin is true, when the chromosome length isn't evenly divided into bin_size-sized bins, 
     if the leftover is less than merge_last_bin_cutoff, merge it into the last bin, otherwise add it as an extra bin 
      (this also ensures that if a chromosome is shorter than a bin, it'll be treated as single bin anyway).
     If normalize_last_bin is also True, normalize the number of positions in the special last bin by its size, 
      so that its value reflects the position density rather than raw count, to match all the other bins for heatmap display.

    If chromosomes is not None, only include those chromosomes; otherwise include all chromosomes in position_dict only.
    """
    # get chromosome list from position_dict - also take care of the case if position_dict keys are (chrom,strand) tuples
    if isinstance(position_dict.keys()[0], str): 
        chromosomes_in_position_dict = position_dict.keys()
    else:                                        
        chromosomes_in_position_dict = set(chrom for (chrom,strand) in position_dict.keys())
    if chromosomes is None:   
        chromosomes = chromosomes_in_position_dict
    # get chromosome lengths, make sure we have them for all chromosomes
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
        chromosomes_no_lengths = set(chromosomes) - set(chromosome_lengths)
        if chromosomes_no_lengths:
            raise Exception("some chromosomes have no length data! %s"%chromosomes_no_lengths)
    chromosome_bin_count_lists = {}
    for chromosome in chromosomes:
        chromosome_length = chromosome_lengths[chromosome]
        # position_dict can be either just chromosome or (chromosome,strand)
        #  - get all the positions on the chromosome in either case, by just adding together the two strand lists if needed.
        #       (convert them to lists first, in case they're numpy arrays, for which + works differently!)
        #       (MAYBE-TODO this list conversion is probably inefficient - could do something with numpy array concatenation...)
        try:                position_list = list(position_dict[chromosome])
        except KeyError:    position_list = list(position_dict[(chromosome,'+')]) + list(position_dict[(chromosome,'-')])
        # divide the chromosome into bin_size-sized ranges, using an x.5 cutoff for clarity
        bin_edges = [x-.5 for x in range(first_bin_offset+1, chromosome_length+1, bin_size)]
        # Make sure there's at least one bin, even if the total length is smaller than a bin!  (for scaffolds/etc)
        # There'll be a smaller-than-bin_size chunk left over at the end (if length doesn't divide evenly into bin_size), 
        #  so add an extra bin for that if it's at least half a bin_size, otherwise make the last bin bigger to compensate.
        if special_last_bin:
            last_bin_edge = chromosome_length-.5
            if len(bin_edges)==1:                               bin_edges.append(last_bin_edge)
            if (last_bin_edge - bin_edges[-1])/bin_size < merge_last_bin_cutoff:  bin_edges[-1] = last_bin_edge
            else:                                                                 bin_edges.append(last_bin_edge)
        # use numpy.histogram to get the actual bin counts, IF we have any bins - otherwise the counts are an empty array
        #  (using numpy.array instead of list just to make sure the types match)
        if len(bin_edges)>1:    bin_count_list, _  = numpy.histogram(position_list, bin_edges)
        else:                   bin_count_list     = numpy.array([])
        # for the last bin in each chromosome, the value should be scaled by the smaller/bigger bin-size...
        if special_last_bin and normalize_last_bin and len(bin_count_list):
            bin_count_list[-1] = bin_count_list[-1] / ((bin_edges[-1] - bin_edges[-2])/bin_size)
        chromosome_bin_count_lists[chromosome] = bin_count_list
    return chromosome_bin_count_lists
    # TODO should probably unit-test this!


def get_20bp_fraction(dataset):
    """ Return the fraction of mutants in dataset that only have 20bp flanking regions (no 21bp ones).

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list of mutants (mutant_analysis_classes.Insertional_mutant instances).
    """
    max_lengths = [max([len(s) for s in mutant.sequences_and_counts.keys()]) for mutant in dataset]
    assert max_lengths.count(20) + max_lengths.count(21) == len(dataset), "Some sequences are outside the 20-21bp range!"
    return max_lengths.count(20)/len(dataset)
    # LATER-TODO should this be here, or a mutant_analysis_classes.Insertional_mutant_pool_dataset method or something?


######### unit-tests

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__(self):
        sys.exit("NO UNIT-TESTS FOR THIS MODULE")
    # LATER-TODO add unit-tests!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
