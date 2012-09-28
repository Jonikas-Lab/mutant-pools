#!/usr/bin/env python
"""
Module containing classes and functions for analysis of deepseq data related to insertional mutant libraries.

This is a module to be imported and used by other programs.  Running it directly runs the built-in unit-test suite.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011
"""

from __future__ import division
# basic libraries
import sys, re
import unittest
from collections import defaultdict
from itertools import combinations
from numpy import median
import copy
# other libraries
import HTSeq
from BCBio import GFF
# my modules
from general_utilities import split_into_N_sets_by_counts, add_dicts_of_ints, keybased_defaultdict, value_and_percentages, FAKE_OUTFILE
from DNA_basic_utilities import SEQ_ENDS, SEQ_STRANDS, SEQ_DIRECTIONS, SEQ_ORIENTATIONS, position_test_contains, position_test_overlap
from seq_basic_utilities import get_seq_count_from_collapsed_header
from deepseq_utilities import check_mutation_count_try_all_methods
# there's a "from parse_annotation_file import parse_gene_annotation_file" in one function, not always needed


class SPECIAL_GENE_CODES(object):
    not_determined = "gene_unknown"
    chromosome_not_in_reference = "unknown_chrom"
    not_found = "no_gene_found"
# it seems like I have to set SPECIAL_GENE_CODES.all_codes afterward because I can't access __dict__ from inside the class
SPECIAL_GENE_CODES.all_codes = [value for (name,value) in SPECIAL_GENE_CODES.__dict__.items() if not name.startswith('__')]

class MutantError(Exception):
    """ Exceptions in the mutant_analysis_classes module; no special behavior."""
    pass


def is_cassette_chromosome(chromosome_name):
    """ Returns True if chromosome_name sounds like one of our insertion cassettes, False otherwise. """
    return ("insertion_cassette" in chromosome_name)

def is_other_chromosome(chromosome_name):
    """ Returns True if chromosome_name is neither cassette nor chromosome/scaffold, False otherwise. """
    if is_cassette_chromosome(chromosome_name):
        return False
    if chromosome_name.startswith('chr') or chromosome_name.startswith('scaffold'):
        return False
    return True

# MAYBE-TODO may want to put the is_*_chromosome functionality under command-line user control someday?  Provide option to give list of cassette chromosome names and "other" ones (or cassette ones and prefixes for genomic ones?)


# MAYBE-TODO it might be good to split this file into multiple files at some point?  At least Insertion_position/etc.

############################ Functions/classes for dealing with alignment/insertion positions ###########################

# has to be a new-style object-based class due to the immutability/hashability thing
class Insertion_position(object):
    """ A descriptor of the position of a genomic insertion, with separate before/after sides; optionally immutable.

    Attributes: 
        - chromosome, strand - the chromosome the insertion is on, and the strand it's in sense orientation to.
        - position_before, position_after - positions before and after the insertion site (1-based): 
                                            integers, or None if unknown.
        - min_position and max_position - lowest/highest possible position values as plain numbers, no ambiguity
        - full_position - string describing the full position: 3-4 for exact positions, 3-? or 4-? if one side is unknown. 
        Note that the three *_position attributes are really property-decorated methods, and cannot be assigned to.

    Methods: 
        - comparison/sorting: __cmp__ is based on chromosome name/number, min_/max_position, strand,
                                   and position_before/_after, in that order - see __cmp__ docstring for more detail
        - copy method returns an identical but separate copy of the object (not just another reference to the same object)
        - printing: __str__ returns a string of the chromosome,strand,full_position values
                    __repr__ returns a string of the Insertion_position() call to create a new identical object.
        Mutability and hashability: by default instances are mutable, and thus unhashable, since they implement __cmp__. 
         There are make_immutable and make_mutable_REMEMBER_CLEANUP_FIRST methods to reversibly toggle the state 
          to immutable/hashable and back. 
         This works by defining __hash__ to use the _make_key() value if immutable and raise an exception otherwise, 
          and decorating __setattr__ and __delattr__ to raise an exception if immutable and work normally otherwise.
         It's not perfect immutability, it can be gotten around by using object.__setitem__ etc, but it's good enough.
    """

    # NOTE: originally I used biopython SeqFeature objects for position_before and position_after (specifically SeqFeature.ExactPosition(min_val) for exactly positions and SeqFeature.WithinPosition(min_val, min_val-max_val) for ambiguous ones), but then I realized I'm not using that and it's over-complicated and hard to make immutable and may not be what I need anyway even if I do need immutable positions, so I switched to just integers. The function to generate those was called make_position_range, and I removed it on 2012-04-23, if I ever want to look it up again later.

    def __init__(self, chromosome, strand, full_position=None, position_before=None, position_after=None, immutable=False):
        """ Initialize all values - chromosome/strand are just copied from arguments; positions are more complicated. 
        
        You must provide either full_position, OR one or both of position_before/position_after. 
        The two position_* arguments must be castable to ints, or None.
        The full_position argument must be a string of the form '100-200', '?-200' or '100-?', such as would be generated 
         by self.full_position() - self.position_before and _after are set based on the two parts of the string.
        Self.min_/max_position are calculated based on self.position_before/_after - both, or whichever one isn't None.
        If immutable is True, the object is made immutable (by calling self.make_immutable() right after initiation.
        """
        # need to make instance mutable to be able to set anything, due to how __setattr__ is decorated
        self.make_mutable_REMEMBER_CLEANUP_FIRST()  
        # now start setting attributes
        self.chromosome = chromosome
        self.strand = strand
        # parse full_position if provided
        if full_position is not None:
            if (position_before is not None) or (position_after is not None):
                raise ValueError("If providing full_position, cannot also provide position_before/position_after!")
            self.position_before, self.position_after = self._parse_full_position(full_position)
        # otherwise use position_before and/or position_after
        else:
            if position_before is None and position_after is None:
                raise ValueError("Can't create an Insertion_position object with no known position values!")
            try:
                self.position_before = None if position_before is None else int(position_before)
                self.position_after = None if position_after is None else int(position_after)
            except TypeError:  
                raise ValueError("position_before/position_after must be int-castable or None!")
        if immutable:   self.make_immutable()

    @property   # this is a builtin decorator to make an attribute out of a method
    def min_position(self):
        if self.position_after is None:     return self.position_before
        elif self.position_before is None:  return self.position_after-1
        else:                               return min(self.position_before, self.position_after-1)

    @property   # this is a builtin decorator to make an attribute out of a method
    def max_position(self):
        if self.position_after is None:     return self.position_before+1
        elif self.position_before is None:  return self.position_after
        else:                               return min(self.position_before+1, self.position_after)

    @property   # this is a builtin decorator to make an attribute out of a method
    def full_position(self):
        info_before = str(self.position_before) if self.position_before is not None else '?'
        info_after = str(self.position_after) if self.position_after is not None else '?'
        return info_before + '-' + info_after

    @classmethod
    def _parse_full_position(cls, full_position_string):
        """ Parse a full_position string to proper (position_before, position_after) value. """
        try:
            before,after = [cls._parse_single_position(s) for s in full_position_string.split('-')]
        except (ValueError,AttributeError):
            raise ValueError("The full_position argument must be a string of the form '100-200', '?-200' or '100-?'!")
        if before is None and after is None:
            raise ValueError("At least one section of the full_position argument must be a number!")
        return before,after

    @staticmethod
    def _parse_single_position(pos_string):
        """ Make a proper position value: cast to int, or return None if '?' or ''. """
        if pos_string in ['?','']:  return None
        else:                       return int(pos_string)

    def __str__(self):
        """ Return short summary of important info. """
        return ' '.join([self.chromosome, self.strand, self.full_position])

    def __repr__(self):
        """ Return full object-construction call (not very complicated) as a string. """
        return "Insertion_position('%s', '%s', full_position='%s', immutable=%s)"%(self.chromosome, self.strand, 
                                                                                   self.full_position, self.immutable)

    def copy(self):
        """ Return a deep-copy of self - NOT just a reference to the same object. """
        return Insertion_position(self.chromosome, self.strand, position_before=self.position_before, 
                                  position_after=self.position_after)

    def _make_key(self):
        """ Make key for sorting/comparison - based on chromosome/position/strand, with improved chromosome-number sorting.

        First two fields are chromosome data - splits chromosome into name/number (both optional), 
         so that "chr2" sorts before "chr12" (but 'chr' before 'chr1', and 'other_chr1' after 'chr4').
        Next two fields are min_/max_position - these are always numerically defined, so ?-101 and 100? will sort together
         (as opposed to if we used position_before/_after, which can be None).
        Next field is strand - we want the sorting on position BEFORE strand, it's more readable/sensible that way.
        Final two fields are position_before/after, to ensure ?-101 isn't considered equal to 100-101.
        """
        chromosome_data = re.search('^(?P<name>.*[^\d])?(?P<number>\d*)', self.chromosome)
        chromosome_name = chromosome_data.group('name')         # if there were no non-digit characters, this is None
        chromosome_number = int(chromosome_data.group('number')) if chromosome_data.group('number') else 0
        all_position_values = (chromosome_name, chromosome_number, self.min_position, self.max_position, 
                               self.strand, self.position_before, self.position_after)
        return all_position_values

    def __cmp__(self,other):
        """ Based on tuple-comparison of (chr_name, chr_number, min_pos, max_pos, strand) - in that order for sane sort."""
        return cmp(self._make_key(), other._make_key())

    # MAYBE-TODO do I also want a rough_comparison method, which would return True for ?-101 and 100-101 etc?  How about 100-101 and 100-102, then?

    # MAYBE-TODO add some kind of merge function to merge two positions into one?  Do I actually need that?  When I'm merging two position-based read groups together because they're actually probably the same mutant with some sequencing indels, I actually just keep the more common position, since that's presumably the correct one. So are there any cases where I'd need to merge positions?

    ### MUTABILITY AND HASHABILITY SWITCHES
    # Sometimes I want to use positions as dictionary keys or put them in sets - so they need to be hashable, 
    #   and since I also need a sane comparison operator, I can't use the default object id-based hashing, 
    #   so I have to make the objects immutable for them to be hashable. 
    # More info on how/why that is: http://docs.python.org/reference/datamodel.html#object.__hash__ 
    # This implementation is not perfectly immutable: you can get around the "immutability" by tricks like in 
    #   make_mutable_REMEMBER_CLEANUP_FIRST, but it's enough to make it clear to the user that it shouldn't be changed.
    # Some links on implementation: http://stackoverflow.com/questions/9997176/immutable-dictionary-only-use-as-a-key-for-another-dictionary, http://stackoverflow.com/questions/1151658/python-hashable-dicts, http://stackoverflow.com/questions/4996815/ways-to-make-a-class-immutable-in-python, http://stackoverflow.com/questions/4828080/how-to-make-an-immutable-object-in-python

    def make_immutable(self):
        """ Reversibly make object immutable (reasonably) and hashable. """
        # just set the flag to make object immutable and hashable
        self.immutable = True

    def make_mutable_REMEMBER_CLEANUP_FIRST(self):
        """ Reversibly make object mutable and non-hashable. REMEMBER TO REMOVE SELF FROM SETS/DICTS BEFORE CALLING! """
        # UNSET the flag to make object immutable and hashable - need to do it in a roundabout way,
        #  because the immutability prevents simply "self.immutable = False" from working!
        self.__dict__['immutable'] = False

    def __hash__(self):
        """ If self.hashable is True, use private _hash method, otherwise raise exception. """
        if  self.immutable:
            return hash(self._make_key())
        else:
            raise MutantError("This %s is currently mutable, and therefore unhashable! "%repr(self)
                              +"Run self.make_immutable() to change this.")

    def exception_if_immutable(function_to_wrap):
        """ Decorator: raise MutantError if self.immutable, else call function as normal. """
        def wrapped_function(self, *args, **kwargs):
            if self.immutable:
                raise MutantError("This %s is currently immutable, cannot change values! "%repr(self)
                                +"You can run self.make_mutable_REMEMBER_CLEANUP_FIRST() to change this - first make SURE "
                                +"to remove it from any sets or dictionary keys that are relying on it being hashable!")
            else:
                return function_to_wrap(self, *args, **kwargs)
        return wrapped_function

    # apply the exception_if_immutable decorator to all methods incompatible with immutability
    __setattr__ = exception_if_immutable(object.__setattr__)
    __delattr__ = exception_if_immutable(object.__delattr__)


def get_insertion_pos_from_HTSeq_read_pos(HTSeq_pos, cassette_end, reads_are_reverse=False, immutable_position=True):
    """ Return a Insertion_position instance giving the cassette insertion position based on HTSeq read position. 

    HTSeq_pos should be a HTSeq.GenomicPosition instance; cassette_end gives the side of the insertion that read is on; 
     reads_are_reverse is True if the read is in reverse orientation to the cassette, False otherwise. 

    The cassette chromosome will be the same as read chromosome; the cassette strand will be the same as read strand, 
     or opposite of the read strand if reads_are_reverse is True.

    The cassette position depends on read position, cassette strand (not read strand) and cassette_end in a complex way:
     Data I have:  which end of the insertion cassette the read is on, and which orientation the cassette is in. 
     Data I want:  the position of the base before and after the insertion, regardless of cassette orientation.
                       (the read orientation in regard to cassette doesn't matter at all here)
     If read is 5' of cassette and cassette is +, or read is 3' of cassette and cassette is -, read is BEFORE cassette
     If read is 3' of cassette and cassette is +, or read is 5' of cassette and cassette is -, read is AFTER cassette
      If read is before cassette, I care about the end of the read; if it's after, I care about the start)
      SAM alignment position is leftmost/rightmost, i.e. end is always the "later" position in the genome, 
       regardless of the read orientation, which gives me what I want, i.e. the insertion position.
    Insertion_position uses a 1-based position system (as opposed to HTSeq, which is 0-based).

    If immutable_position is True, the position will be made immutable after creation (this is reversible).
    """
    ### chromosome is always the same as read (also checking that the HTSeq_pos has a chrom attribute at all here)
    try:
        chrom = HTSeq_pos.chrom
    except AttributeError:
        raise ValueError("Invalid position %s! Need an HTSeq iv object. (If empty, maybe read wasn't aligned?)"%HTSeq_pos)
    ### cassette strand is the same as read strand, OR the opposite if reads_are_reverse is True
    strand = HTSeq_pos.strand
    if reads_are_reverse:     strand = '+' if strand=='-' else '-'
    ### cassette position depends on the read position and cassette_end in a somewhat complex way 
    #   (see description in docstring, and ../notes.txt for even more detail)
    # HTSeq is 0-based and I want 1-based, thus the +1; end has no +1 because in HTSeq end is the base AFTER the alignment.
    if (cassette_end=='5prime' and strand=='+') or (cassette_end=='3prime' and strand=='-'):      
        pos_before, pos_after = HTSeq_pos.end, None
    elif (cassette_end=='3prime' and strand=='+') or (cassette_end=='5prime' and strand=='-'):    
        pos_before, pos_after = None, HTSeq_pos.start+1
    else:                           
        raise ValueError("cassette_end argument must be one of %s."%SEQ_ENDS)
    return Insertion_position(chrom, strand, position_before=pos_before, position_after=pos_after, 
                              immutable=immutable_position)


def find_gene_by_pos(insertion_pos, chromosome_GFF_record, detailed_features=False, quiet=False):
    """ Look up insertion_pos in chromosome_GFF_record; return (gene_ID,orientation,subfeature) of insertion in gene.

    Insertion_pos is an Insertion_position instance.  If insertion_pos overlaps a gene in chromosome_GFF_record, 
     return geneID (with '(edge)' appended if insertion_pos isn't completely contained inside the gene), 
      orientation ('sense' if insertion_pos and gene are in the same direction, 'antisense' otherwise), 
      and the name of the subfeature (exon/intron/UTR) Insertion_pos is in (or '?' if detailed_features is False); 
      if insertion_pos is not in a gene, return ('no_gene_found', '-', '-').
     If Insertion_pos is on the edge of two features, subfeature will be 'X/Y'; if it's something more 
      unlikely/complicated, subfeature will be 'X/Y/Z??' and a warning will be printed to STDOUT.

    Chromosome_GFF_record is a record generated by BCBio.GFF parser from a gff file (usually one file = multiple records). 
    """
    # see notes_on_GFF_parsing.txt for what a GFF record (chromosome_GFF_record) will be like
    assert insertion_pos.strand in ['+','-','both'], "Strand should be %s, and is %s!"%(' or '.join(SEQ_STRANDS),
                                                                                       insertion_pos.strand)
    # go over all the genes in the chromosome record and look for one that matches the position in insertion_pos
    ins_start,ins_end = insertion_pos.min_position, insertion_pos.max_position
    for gene in chromosome_GFF_record.features:
        # for GFF positions, always add 1 to the gene/feature start, because BCBio uses 0-based and I use 1-based, 
        #  but don't add 1 to the end, because BCBio uses end-exclusive and I use end-inclusive.
        gene_start, gene_end = gene.location.start.position+1, gene.location.end.position
        if position_test_overlap(gene_start, gene_end, ins_start,ins_end):
            # if insertion_pos is inside the gene, save the gene ID
            gene_ID = gene.id
            # if it overlaps an edge, note that in features_gene by adding 'gene_edge'
            # (MAYBE-TODO if I look at things like inner/outer flanking regions, this is where that would go as well)
            if position_test_contains(gene_start, gene_end, ins_start,ins_end): features_gene = []
            else:                                                               features_gene = ['gene_edge']
            # calculate orientation of insertion_pos vs gene
            if insertion_pos.strand=='both':    orientation = 'both'
            elif gene.strand==1:                orientation = 'sense' if insertion_pos.strand=='+' else 'antisense'
            elif gene.strand==-1:               orientation = 'sense' if insertion_pos.strand=='-' else 'antisense'
            else:                               orientation = '?'
            # if we're not looking for detailed features, just return the data now, using '?' as gene_feature
            if not detailed_features:
                full_feature = '/'.join(features_gene+['?'])
                return gene_ID, orientation, full_feature

            ### Find which feature of the gene the insertion_pos should be annotated as:
            # if gene has no features listed, use 'no_mRNA' as gene_feature (different from '?' or '-')
            if len(gene.sub_features)==0:
                inner_feature = 'no_mRNA'
            # if gene feature/subfeature structure isn't as expected, print warning and return '??'
            elif len(gene.sub_features)>1:
                if not quiet:
                    print("Warning: gene %s in gff file has multiple sub-features (mRNAs?)! "%gene_ID
                          +"Returning '??' feature.")
                inner_feature = '??'
            elif gene.sub_features[0].type != 'mRNA':
                if not quiet:
                    print("Warning: gene %s in gff file has unexpected non-mRNA sub-features! "%gene_ID
                          +"Returning '??' feature.")
                inner_feature = '??'
            else:
                mRNA = gene.sub_features[0]
                mRNA_start, mRNA_end = mRNA.location.start.position+1,mRNA.location.end.position
                # if insertion_pos is outside the mRNA, use 'outside_mRNA' as inner_feature
                if not position_test_overlap(mRNA_start, mRNA_end, ins_start, ins_end):
                    inner_feature = 'outside_mRNA'
                # if insertion_pos is inside the mRNA and mRNA has no subfeatures, use 'mRNA_no_exons' as inner_feature
                elif len(mRNA.sub_features)==0:   
                    if position_test_contains(mRNA_start, mRNA_end, ins_start, ins_end):
                        inner_feature = 'mRNA_no_exons'
                    else:
                        inner_feature = 'mRNA_edge'
                else: 
                    # otherwise go over all subfeatures, see which ones contain/overlap insertion_pos
                    #  (check for overlap only if the contains test failed)
                    features_inside = []
                    if position_test_contains(mRNA_start, mRNA_end, ins_start, ins_end): features_edge = []
                    else:                                                                features_edge = ['mRNA_edge']
                    for feature in mRNA.sub_features:
                        feature_start, feature_end = feature.location.start.position+1, feature.location.end.position
                        if position_test_contains(feature_start, feature_end, ins_start, ins_end):
                            features_inside.append(feature.type)
                        elif position_test_overlap(feature_start, feature_end, ins_start, ins_end):
                            features_edge.append(feature.type)
                    # MAYBE-TODO may want to treat exons before 5'UTR or after 3'UTR specially? (EH, none in current file)
                    # MAYBE-TODO may want to treat cases with multiple UTRs specially?  We do have those in current file!
                    # if insertion_pos is inside a single mRNA subfeature, use the type of the subfeature as inner_feature
                    if len(features_inside)==1 and len(features_edge)==0:
                        inner_feature = features_inside[0]
                    # if insertion_pos is on the edge of two mRNA subfeatures, use 'subfeature1/subfeature2'
                    elif len(features_inside)==0 and len(features_edge)==2:
                        inner_feature = '/'.join(features_edge)
                        # MAYBE-TODO treat insertions CLOSE to an edge specially too? How large is a splice junction?
                    # if insertion_pos is on the edge of ONE mRNA subfeature, or not touching any subfeatures at all, 
                    #  the implied subfeature is either an intron (if between features) or mRNA_before/after_exons, 
                    #   (which shouldn't happen in normal files).
                    elif len(features_inside)==0 and len(features_edge)<=1:
                        # figure out what the implied feature is
                        if ins_start < min([feature.location.start.position+1 for feature in mRNA.sub_features]):
                            implied_feature = 'mRNA_before_exons'
                        elif ins_end > max([feature.location.end.position for feature in mRNA.sub_features]):
                            implied_feature = 'mRNA_after_exons'
                        else:
                            implied_feature = 'intron'
                        # set inner_feature based on whether insertion_pos is on a real/implied feature edge 
                        #  or completely inside an implied feature
                        if len(features_edge)==1:
                            inner_feature = features_edge[0] + '/' + implied_feature
                        elif len(features_edge)==0:
                            inner_feature = implied_feature
                    # if insertion_pos is inside two features, or inside one and on the edge of another, 
                    #  print a warning, and use all the feature names, with a ?? at the end to mark strangeness
                    else:
                        inner_feature = '/'.join(features_inside+features_edge) + '??'
                        if not quiet:
                            print(("Warning: Location (%s,%s) matched multiple features (%s) "
                                  +"in gene %s!")%(ins_start, ins_end, inner_feature, gene_ID)) 
                # MAYBE-TODO also output distance from start/end of gene/feature?
            
            # prepend whatever gene-level features (edge etc, or []) were found at the start to the full value
            full_feature = '/'.join(features_gene+[inner_feature])
            return gene_ID, orientation, full_feature
    # MAYBE-TODO do I want to consider the case of an insertion on the edge between two genes?  Or even in two genes and not on the edge, if there are overlapping genes, but hopefully there aren't!  (No overlapping or immediately adjacent genes found in our current file - lowest between-gene distance is 1)
    # if no gene matching insertion_pos was found, return special value
    return SPECIAL_GENE_CODES.not_found, '-', '-'
    # MAYBE-TODO add unit tests?  But this is included in a pretty thorough run-test, so may not be necessary.


############################### Main classes describing the mutants and mutant sets #####################################

# help functions returning "blank" mutants for defaultdicts (can't use lambdas because pickle doesn't like them)
def blank_readcount_only_mutant():
    return Insertional_mutant(readcount_related_only=True)
def blank_full_single_mutant():
    return Insertional_mutant(multi_dataset=False,readcount_related_only=False)
def blank_single_mutant_with_pos(pos):
    return Insertional_mutant(insertion_position=pos, multi_dataset=False, readcount_related_only=False)
def blank_multi_mutant_with_pos(pos):
    return Insertional_mutant(insertion_position=pos, multi_dataset=True, readcount_related_only=False)


class Insertional_mutant():
    """ Data regarding a particular insertional mutant: insertion position, gene, read numbers/sequences, etc.

    Mutants have the following attributes (data):
     1) Position/gene attributes (except readcount-related-only mutants, which don't have those):
       - position - an Insertion_position instance giving the insertion chromosome/strand/position
       - gene, orientation, gene feature - what gene the insertion is in (or one of the SPECIAL_GENE_CODES if unknown), 
                   whether it's in the sense or antisense orientation vs the gene, what feature (exon/intron/UTR) it's in.
     2) Readcount-related attributes (multi-dataset mutants don't have those on the top-level):
       - total_read_count, perfect_read_count - number of all and perfectly aligned deepseq reads
       - unique_sequence_count, sequences_and_counts - number of unique read sequences, and a seq:count dictionary
     3) Attributes related to multi-dataset mutants:
       - multi_dataset - a True/False flag specifying whether a mutant is multi-dataset or not
       - by_dataset - a dataset_name:readcount_data_object dictionary which multi-dataset mutants contain 
           instead of top-level readcount-related attributes. See "Types of mutants" section for more info/examples.

    Types of mutants (with corresponding flags to the __init__ function):
       - normal or readcount-related-only - readcount-related-only mutants have no position/gene attributes, 
           only readcount-related attributes. Used as parts of a multi-dataset mutants to store per-dataset data.
       - single-dataset or multi-dataset - the multi_dataset attribute is True if multi-dataset, false otherwise.
         Both versions have one copy of position/gene attributes. A single-dataset mutant also has just one copy of 
          all readcount-related attributes (and doesn't store the dataset name for them); a multi-dataset mutant has 
          a by_dataset dataset_name:readcount_only_mutant dictionary instead of top-level readcount-related attributes; 
           the dictionary contains a separate readcount-related-only mutant for each dataset, 
            which stores all the readcount-related attributes for that dataset.
          For example to get total_read_count from a single-dataset mutant, you use mutant.total_read_count; 
           to get it for dataset d1 from a multi-dataset mutant, you use mutant.by_dataset['d1'].total_read_count.
         A single-dataset mutant can be converted to a multi-dataset one in place using ____; 
         a new single-dataset mutant can be extracted from a multi-dataset one using _____. 

    Methods (functions) of mutant objects:
        For detailed information on mutant methods, see method docstrings.
        Methods with names starting with _ are private and shouldn't be used from outside the object itself.
        Readcount-related methods take a dataset_name argument - when calling those methods on a single-dataset mutant, 
         don't give a dataset_name value, but when calling them on multi-dataset mutant, give the name of the dataset you 
          want to apply the method to.  Doing this wrong will raise an MutantError.  
         Some methods (like get_main_sequence) have a well-defined behavior when called on a multi-dataset mutant without 
          specifying a dataset - they'll give the result for all datasets added together.
         Some methods may not have multi-dataset functionality implemented, if I didn't think it would be useful.
    """
    # MAYBE-TODO rewrite this with subclassing, so that multi-dataset mutants don't HAVE some methods (only the read-related sub-mutants have these methods, etc?  Instead of adding a dataset arg to each method to make it work on both multi and single, just require the caller to call mutant.read_info(dataset).method for multi-mutants, for read-related methods.

    def __init__(self, insertion_position=None, multi_dataset=False, readcount_related_only=False):
        """ Set self.position based on argument; initialize read/sequence counts to 0 and gene-info to unknown. 

        insertion_position argument should be a Insertion_position instance. If multi_dataset is True, all 
         readcount-related attributes will be initialized as dataset_name:value defaultdicts instead of single values.
        If readcount_related_only is True, non-readcount-related attributes won't be initialized at all.
        """
        self.multi_dataset = multi_dataset
        # "standard" (not readcount-related-only) mutants need general attributes:
        if not readcount_related_only:
            self.position = insertion_position
            # MAYBE-TODO should I have a class for the gene data? Especially if I add more of it (annotation info etc)
            self.gene = SPECIAL_GENE_CODES.not_determined
            self.orientation = '?'
            self.gene_feature = '?'
        # single-dataset mutants get readcount-related attributes; multi-dataset mutants get a by_dataset dictionary, 
        #  with dataset names as keys and readcount-related-only mutants as values.
        if not multi_dataset:
            self._set_readcount_related_data_to_zero()
        else:
            self.by_dataset = defaultdict(blank_readcount_only_mutant)

    def read_info(self, dataset=None, strict=False):
        """ Help function to get read-info-containing object for both multi-dataset and single (dataset=None) mutants.

        For single-dataset (dataset=None), return self.
        For multi-dataset, return self.by_dataset[dataset] if present - if not present, raises an exception if strict, 
         otherwise returns an empty read-info object.
        """
        if dataset is None and self.multi_dataset:  
            raise MutantError("This is a multi-dataset mutant - must provide dataset arg!")
        if dataset is not None and not self.multi_dataset:  
            raise MutantError("This is not a multi-dataset mutant - cannot provide dataset arg!")
        if self.multi_dataset and dataset not in self.by_dataset and strict:
            raise MutantError("No dataset %s in this multi-dataset mutant! Present datasets are %s"%(dataset, 
                                                                                                 self.by_dataset.keys()))
        if self.multi_dataset:  
            try:                return self.by_dataset[dataset]
            except KeyError:    return blank_readcount_only_mutant()
        else:                   return self
        # TODO unit-tests?

    # MAYBE-TODO give each mutant some kind of unique ID at some point in the process?  Or is genomic location sufficient?  If we end up using per-mutant barcodes (in addition to the flanking sequences), we could use that, probably, or that plus genomic location.

    def _set_readcount_related_data_to_zero(self):
        """ Set all readcount-related data to 0/empty."""
        if self.multi_dataset:
            raise MutantError("This is a multi-dataset mutant - no top-level readcount data!")
        self.total_read_count      = 0
        self.perfect_read_count    = 0
        self.unique_sequence_count = 0
        self.sequences_and_counts  = defaultdict(int)

    def _check_consistent_multi_dataset_args(self, dataset_name, function_name='<some_func>', multi_noname_allowed=False):
        """ Make sure dataset_name is consistent with single-/multi-dataset status of mutant; raise MutantError if not.

        Raise MutantError if dataset_name is None for multi-dataset mutant (i.e. the caller is trying to do some operation 
          on a multi-dataset mutant without specifying which dataset to do it on (there are a few cases in which this
           does make sense - in that case multi_noname_allowed=True should be passed), 
          or if dataset_name is not None for single-dataset mutant (i.e. the caller is trying to do some operation 
           on a specific dataset when the mutant is single-dataset and has no named datasets at all).
        The exception text will use function_name as the name of the calling function, if given.
        """
        # MAYBE-TODO probably possible to do this as a decorator and without passing function_name explicitly, but I'm not sure how...  See http://stackoverflow.com/questions/5063607/is-there-a-self-flag-can-reference-python-function-inside-itself for how to get the function name, but I don't actually know how to do the rest with a decorator (specifically how to get the dataset_name value from the decorated function so I can do things with it!)
        if self.multi_dataset and (dataset_name is None) and not multi_noname_allowed: 
            raise MutantError("This mutant is in multi-dataset form! Provide dataset_name "
                              +"to specify which dataset to apply %s() to."%function_name)
        if not self.multi_dataset and (dataset_name is not None): 
            raise MutantError("You're trying to apply %s() to dataset %s, but this mutant "%(function_name, dataset_name)
                              +"hasn't been converted to multi-dataset form! Run convert_to_multi_dataset first.")

    def add_read(self, HTSeq_alignment, read_count=1, treat_unknown_as_match=False, dataset_name=None):
        """ Add a read to the data (or multiple identical reads, if read_count>1); return True if perfect match.

        Specifically: increment total_read_count, increment perfect_read_count if read is a perfect 
         alignment, increment the appropriate field of sequences_and_counts based on read sequence.
        Note: this does NOT check the read position to make sure it matches that of the object.

        If self.multi_dataset is True, dataset_name is required and the appropriate by-dataset values will be changed.
        """
        ### standard stuff for multi-dataset mutants
        # make sure the caller isn't confused about single/multi-dataset input
        self._check_consistent_multi_dataset_args(dataset_name, 'add_read')
        # if it's a single-dataset mutant, all the readcount-related data is in self directly;
        #  if it's a multi-dataset mutant, the readcount-related data we want is in self.by_dataset[dataset_name]
        if not self.multi_dataset:  readcount_data_container = self
        else:                       readcount_data_container = self.by_dataset[dataset_name]
        # MAYBE-TODO these first three lines show up in a lot of functions here - could I make them a decorator?

        ### main functionality
        # MAYBE-TODO check HTSeq_alignment chromosome/strand to make sure it matches data in self?  Don't check position, that's more complicated (it can be either start or end) - could maybe check that position is within, idk, 10bp of either alignment start or alignment end?  Or not - I may want to cluster things in a non-position-based way anyway!  Hmmm...
        seq = HTSeq_alignment.read.seq
        # if it's a new sequence, increment unique_sequence_count; add a count to the sequences_and_counts dictionary.
        if seq not in readcount_data_container.sequences_and_counts:
            readcount_data_container.unique_sequence_count += 1
        readcount_data_container.sequences_and_counts[seq] += read_count
        # increment total_read_count
        readcount_data_container.total_read_count += read_count
        # figure out if the read is perfect and increment perfect_read_count if yes; return True if perfect else False.
        treat_unknown_as = 'match' if treat_unknown_as_match else 'mutation'
        mutation_count = check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as=treat_unknown_as)
        if mutation_count==0:  
            readcount_data_container.perfect_read_count += read_count
            return True
        else:
            return False

    def update_gene_info(self, gene, orientation, gene_feature):
        """ Update gene/orientation/feature: if both are the same or one is unknown, keep known; if different, raise error.
        """
        # grab the set of non-special values from own and new data
        gene_both = set([self.gene, gene]) - set([SPECIAL_GENE_CODES.not_determined])
        orientation_both = set([self.orientation, orientation]) - set(['?'])
        feature_both = set([self.gene_feature, gene_feature]) - set(['?'])
        # if there are two different non-special values for any of the three attributes, raise MutantError
        if not all([len(data_both)<=1 for data_both in [gene_both, orientation_both, feature_both]]):
            raise MutantError("Can't merge the two mutants: the gene/orientation/feature data differs!")
        # otherwise set own data to the better one of own/new data (unless neither are present)
        if gene_both:           self.gene = gene_both.pop()
        if orientation_both:    self.orientation = orientation_both.pop()
        if feature_both:        self.gene_feature = feature_both.pop()

    def merge_mutant(self, other, check_gene_data=True, opposite_strand_tandem=False):
        """ Merge other mutant into this mutant: merge counts, sequences, optionally position; set other's counts to 0.
        Does NOT check that the positions match - this should be done by the caller. 
        """
        if self.multi_dataset or other.multi_dataset:  
            raise MutantError("Merging multi-dataset mutants not implemented!")
            # MAYBE-TODO implement merging for multi-dataset mutants?  Seems like unnecessary complication for now.

        # make sure the two mutants don't have conflicting gene data, if required  (and update if own gene data is unknown)
        # (note: NOT checking position) (this needs to be done first, BEFORE we make changes to the mutants!)
        if check_gene_data:
            try:
                # in the opposite_strand_tandem case, we expect orientations to be different, so pass '?' to the check
                if opposite_strand_tandem:  self.update_gene_info(other.gene, '?', other.gene_feature)
                else:                       self.update_gene_info(other.gene, other.orientation, other.gene_feature)
            except MutantError:
                raise MutantError("Can't merge the two mutants: the gene/orientation/feature data differs!")
        
        ### merge positions
        # in the opposite_strand_tandem case, the new strand should be "both" and the new position should give both ends
        if opposite_strand_tandem:
            assert self.position.strand!=other.position.strand
            # keep track of original by-strand readcounts
            self.original_strand_readcounts = dict([(m.position.strand, m.total_read_count) for m in (self, other)])
            self.position.strand = 'both'
            self.position.position_before = self.position.position_before or other.position.position_before
            self.position.position_after = self.position.position_after or other.position.position_after
            # also set the gene-orientation to "both", unless there is no gene
            if self.gene not in SPECIAL_GENE_CODES.all_codes:
                self.orientation = 'both'
        # otherwise, just keep the self position, since that should be the one with more reads; 
        else:
            assert self.total_read_count >= other.total_read_count, "Merging mutants the wrong way!"
            # MAYBE-TODO implement position merging?  Are we ever going to need it, really?  PROBABLY NOT.

        ### merge read counts
        self.total_read_count += other.total_read_count
        # increment perfect read counts only for opposite-side-tandem, because in that case both sides are equally "good"; 
        #  in cases of merging due to off-by-one position, "perfect" read counts from the wrong position aren't perfect!
        if opposite_strand_tandem:
            self.perfect_read_count += other.perfect_read_count
        # merge sequences
        for (seq,count) in other.sequences_and_counts.iteritems():
            self.sequences_and_counts[seq] += count
        self.unique_sequence_count = len(self.sequences_and_counts)
        # LATER-TODO may want to keep more data about sequences! Like exact position and strand and number of mutations - may want to store a list of HTSeq.alignment objects instead of just sequences+counts, really
        # TODO for merged opposite-tandems, should keep info about merging in the mutant itself!  How, exactly?  I did already add m.original_strand_readcounts...
        # MAYBE-TODO for merged adjacent mutants, keep track of merging info too?
        other._set_readcount_related_data_to_zero()

    # MAYBE-TODO should there also be an add_mutant function that adds mutant readcounts together and optionally makes sure the positions are the same, or should that be the same as merge_mutant?  There are two separate use cases: one where we're merging two adjacent mutants from one dataset, one where we're adding the counts for two mutants from two different datasets. But the mechanics are similar...

    def add_counts(self, total_count, perfect_count, sequence_variant_count, assume_new_sequences=False, 
                   dataset_name=None):
        """ Increment self.total_read_count, self.perfect_read_count and self.unique_sequence_count based on inputs.

        Note that if self.unique_sequence_count>0, it's impossible to determine the correct new value: 
         if we had old data with one unique sequence and now we have new data with another one, how do we know
          if that's the same or different sequence?  The correct total could be 1 or 2, so it's an option:
         If assume_new_sequences is True, the total is old+new; if it's False, the total is max(old,new).

        If self.multi_dataset is True, dataset_name is required and the appropriate by-dataset values will be changed.
        """
        self._check_consistent_multi_dataset_args(dataset_name, 'add_counts')
        if not self.multi_dataset:  readcount_data_container = self
        else:                       readcount_data_container = self.by_dataset[dataset_name]

        readcount_data_container.total_read_count += total_count
        readcount_data_container.perfect_read_count += perfect_count
        if assume_new_sequences:
            readcount_data_container.unique_sequence_count += sequence_variant_count
        else:
            readcount_data_container.unique_sequence_count = max(readcount_data_container.unique_sequence_count,
                                                                 sequence_variant_count)

    def add_sequence_and_counts(self, seq, seq_count, add_to_uniqseqcount=True, dataset_name=None):
        """ Add seq_count to self.sequences_and_counts[seq] (it's created with count 0 if seq wasn't a key before).

        Note: if add_to_uniqseqcount is False, this will never increment self.unique_sequence_count;
         otherwise it only does so if seq was not already present in the self.sequences_and_counts data
          and if the total number of sequences in self.sequences_and_counts is higher than self.unique_sequence_count. 

        If self.multi_dataset is True, dataset_name is required and the appropriate by-dataset values will be changed.
        """
        self._check_consistent_multi_dataset_args(dataset_name, 'add_sequence_and_counts')
        if not self.multi_dataset:  readcount_data_container = self
        else:                       readcount_data_container = self.by_dataset[dataset_name]
        # increment unique_sequence_count if desired (needs to be done first because of the checks it's doing)
        if add_to_uniqseqcount:
            if seq not in readcount_data_container.sequences_and_counts and\
               len(readcount_data_container.sequences_and_counts)>readcount_data_container.unique_sequence_count:
                    readcount_data_container.unique_sequence_count += 1
        # main function: add another seq_count counts of seq
        readcount_data_container.sequences_and_counts[seq] += seq_count

    def get_main_sequence(self, N=1, dataset_name=None):
        """ Return the most common sequence in this mutant and its count (or Nth most common sequence if N is provided).

        If mutant is multi-dataset and dataset_name is given, return the most common sequence for just that dataset; 
         or if dataset_name is None, return most common sequence by total count over all the datasets.
        """
        # MutantError if dataset_name given for single-dataset mutant, but not giving it for multi-dataset is allowed
        self._check_consistent_multi_dataset_args(dataset_name, 'get_main_sequence', multi_noname_allowed=True)
        ### pick the seq:count dictionary to use:
        # if the mutant is single-dataset, just return top sequence
        if not self.multi_dataset:
            seqs_to_counts = self.sequences_and_counts
        # if mutant is multi-datset, if dataset_name was given, just take that, otherwise take sum of all datasets
        elif self.multi_dataset:
            if dataset_name is not None:
                seqs_to_counts = self.by_dataset[dataset_name].sequences_and_counts
            else:
                seqs_to_counts = reduce(add_dicts_of_ints, 
                                        [dataset.sequences_and_counts for dataset in self.by_dataset.itervalues()])
                # MAYBE-TODO print a warning if different dataset mutants have different main sequences?

        ### now actually find the Nth main sequence in the chosen seq:count dictionary
        sequences_by_count = sorted([(count,seq) for (seq,count) in seqs_to_counts.iteritems()], reverse=True)
        # try returning the Nth sequence and count; return nothing if there are under N sequences.
        try:                return tuple(reversed(sequences_by_count[N-1]))
        except IndexError:  return ('',0)
        # MAYBE-TODO should probably make that '-' or something instead of '', empty strings are hard to see. 
        #  On the other hand '-' isn't a valid sequence, and '' is...
        # MAYBE-TODO should there be a warning/failure/something if it's a multi-dataset mutant and the user wants
        #  an overall main sequence and only some of the mutants have any sequence data?

    def _copy_non_readcount_data(self, source_mutant):
        """ Copy non-readcount-related data from source_mutant to self (making new copies of all objects). """
        # COPY the position, not just make another name for the same value - I wrote a copy() function for positions
        self.position     = source_mutant.position.copy() 
        # strings are immutable and thus safe to "copy" by adding another name to the same value
        self.gene         = source_mutant.gene
        self.orientation  = source_mutant.orientation
        self.gene_feature = source_mutant.gene_feature

    def _copy_readcount_related_data(self, source_mutant):
        """ Copy readcount-related data from source_mutant to self (making new copies of all objects). """
        # integers are immutable and thus safe to "copy" by adding another name to the same value
        self.total_read_count      = source_mutant.total_read_count
        self.perfect_read_count    = source_mutant.perfect_read_count
        self.unique_sequence_count = source_mutant.unique_sequence_count
        # using dict to make a COPY of the dict instead of just creating another name for the same value
        self.sequences_and_counts  = dict(source_mutant.sequences_and_counts)

    def convert_to_multi_dataset(self, current_dataset_name=None, ignore_if_already_multi=False):
        """ Convert mutant from single-dataset to multi-dataset mutant; assign readcount-data to current_dataset_name.
        
        Set self.multi_dataset to True, set up the by_dataset dictionary; remove all current readcount-related data
         (move it to self.by_dataset dictionary[current_dataset_name] first if current_dataset_name is not None).
        The position and gene information will remain single, since it shouldn't vary between datasets.
        """
        if self.multi_dataset:
            if ignore_if_already_multi:
                return
            else:
                raise MutantError("Can't run convert_to_multi_dataset - mutant is already multi-dataset!")
        # generate new by-dataset dictionary of all the readcount-related data separately for each dataset
        self.by_dataset = defaultdict(blank_readcount_only_mutant)
        self.multi_dataset = True
        # if current_dataset_name isn't None, copy the old readcount-related data values from self 
        #  to the new dictionaries using the add_other_mutant_as_dataset with self as other_mutant
        if current_dataset_name is not None:
            self.add_other_mutant_as_dataset(self, current_dataset_name, check_constant_data=False)
        # remove all old single-dataset readcount-related data
        del self.total_read_count
        del self.perfect_read_count
        del self.unique_sequence_count
        del self.sequences_and_counts

    def add_other_mutant_as_dataset(self, other_mutant, other_mutant_dataset_name, 
                                    overwrite=False, check_constant_data=False):
        """ Copy all readcount-related data from other_mutant to self.by_dataset dictionary[other_mutant_dataset_name].

        If self isn't a multi-dataset mutant, raise an exception.
        If check_constant_data is True, check that the position/gene data of self and other_mutant matches.
        If self already has a other_mutant_dataset_name dataset, raise MutantError, unless overwrite=True, then overwrite.
        """
        if not self.multi_dataset:
            raise MutantError("This mutant hasn't been converted to multi-dataset form, "
                              +"can't run add_other_mutant_as_dataset! Run convert_to_multi_dataset first.")
        if other_mutant_dataset_name in self.by_dataset and not overwrite:
            raise MutantError("This mutant already has a %s dataset! Can't overwrite it with "%other_mutant_dataset_name
                              +"new one.  Choose a different name for new dataset, or use overwrite=True argument.")

        # if desired, check that the position/gene data matches (and update if own gene data is unknown)
        #  (probably should be using ifs rather than asserts, but I think since they're wrapped in a try/except it's fine)
        if check_constant_data:
            try:
                # MAYBE-TODO should this be an inexact comparison? Right now 100-101 is NOT equal to ?-101 or to 100-102.
                assert self.position == other_mutant.position
                self.update_gene_info(other_mutant.gene, other_mutant.orientation, other_mutant.gene_feature)
            except (AssertionError,MutantError):
                raise MutantError("Can't add mutant2 as dataset to mutant1: the mutant position/gene data differs!")

        # make a new empty Insertional_mutant object to hold the readcount-related data from other_mutant, 
        #  and put it in the self.by_dataset dictionary under other_mutant_dataset_name
        self.by_dataset[other_mutant_dataset_name] = Insertional_mutant(readcount_related_only=True)
        # now fill this new object with readcount-related data from other_mutant
        self.by_dataset[other_mutant_dataset_name]._copy_readcount_related_data(other_mutant)

    def give_single_dataset_mutant(self, single_dataset_name, force=False):
        """ Return a single-dataset mutant based on single_dataset_name; don't modify current mutant.

        If there is no single_dataset_name in current mutant's by_dataset dictionary, raise exception, 
         unless force is True, then return new mutant with zero read-count.
        """
        if not self.multi_dataset:
            raise MutantError("This mutant is not multi-dataset, can't run give_single_dataset_mutant!")
        if single_dataset_name not in self.by_dataset.keys() and not force:
            raise MutantError("This mutant doesn't have a %s dataset! "%single_dataset_name
                              +"Use force=True argument if you want a zero-readcount mutant returned anyway.")
        # generate new mutant, fill it with readcount-related data from self.by_dataset[single_dataset_name] 
        #  and general data from self
        new_mutant = Insertional_mutant(multi_dataset=False)
        new_mutant._copy_non_readcount_data(self)
        new_mutant._copy_readcount_related_data(self.by_dataset[single_dataset_name])
        return new_mutant

    def give_all_single_dataset_mutants(self):
        """ Split multi-dataset mutant into a dataset_name:single-dataset_mutant dictionary and return it; 
        don't modify original mutant.
        """
        if not self.multi_dataset:
            raise MutantError("This mutant is not multi-dataset, can't run give_all_single_dataset_mutants!")
        mutant_dictionary = defaultdict(blank_full_single_mutant)
        for dataset_name in self.by_dataset.iterkeys():
            mutant_dictionary[dataset_name] = self.give_single_dataset_mutant(dataset_name)
        return mutant_dictionary

    # MAYBE-TODO should there be a copy_mutant function to make a deepcopy?


class Dataset_summary_data():
    """ Summary data for a Insertional_mutant_pool_dataset object.  Lots of obvious attributes; no non-default methods.
    """
    # LATER-TODO update docstring

    def __init__(self, dataset, cassette_end, reads_are_reverse, dataset_name=None):
        """ Initialize everything to 0/empty/unknown. """
         # make sure the arguments are valid values
        if not cassette_end in SEQ_ENDS+['?']: 
            raise ValueError("The cassette_end variable must be one of %s or '?'!"%SEQ_ENDS)
        if not reads_are_reverse in [True,False,'?']: 
            raise ValueError("The reads_are_reverse variable must be True, False, or '?'!")
        # reference to the containing dataset (for read-counting purposes etc), 
        #  and the dataset name (None if it's a single dataset, string for multi-datasets)
        self.dataset_name = dataset_name
        self.dataset = dataset
        # information on reads that aren't included in the dataset mutants
        self.discarded_read_count, self.discarded_wrong_start, self.discarded_no_cassette = 'unknown', 'unknown', 'unknown'
        self.non_aligned_read_count, self.unaligned, self.multiple_aligned = 0, 'unknown', 'unknown'
        self.ignored_region_read_counts = defaultdict(int)
        # mutant merging information
        NaN = float('NaN')
        self.adjacent_max_distance = None
        self.adjacent_merging_count_ratio = None
        self.merging_which_chromosomes = (None, None)
        self.merged_adjacent_pairs = 0
        self.merged_opposite_tandems = 0
        # TODO should the adjacent-mutant-counts (same_position_opposite, adjacent_opposite_toward, adjacent_opposite_away, adjacent_same_strand) be constants or property values?  May want to keep them as constants, they're expensive to calculate...
        self.same_position_opposite = NaN
        self.adjacent_opposite_toward = NaN
        self.adjacent_opposite_away = NaN
        self.adjacent_same_strand = NaN
        # MAYBE-TODO should cassette_end and reads_are_reverse be specified for the whole dataset, or just for each set of data added, in add_alignment_reader_to_data? The only real issue with this would be that then I wouldn't be able to print this information in the summary - or I'd have to keep track of what the value was for each alignment reader added and print that in the summary if it's a single value, or 'varied' if it's different values. Might also want to keep track of how many alignment readers were involved, and print THAT in the summary!  Or even print each (infile_name, cassette_end, reads_are_reverse) tuple as a separate line in the header.
        self.cassette_end = cassette_end
        self.reads_are_reverse = reads_are_reverse

    # TODO unit-test all these methods!

    def add_discarded_reads(self, N_all_discarded, N_wrong_start, N_no_cassette, replace=False):
        """ Add not-None arg values to discarded_read_count, discarded_wrong_start and discarded_no_cassette (or replace). 
        
        If the original values are 'unknown', or replace is True, replace instead of adding.
        If any of the args is None, don't modify the original value, unles replace is True, then set to 'unknown'.
        """
        if N_all_discarded is not None:
            if replace or self.discarded_read_count=='unknown': self.discarded_read_count = int(N_all_discarded)
            else:                                               self.discarded_read_count += int(N_all_discarded)
        elif replace:                                           self.discarded_read_count = 'unknown'
        if N_wrong_start is not None:
            if replace or self.discarded_wrong_start=='unknown': self.discarded_wrong_start = int(N_wrong_start)
            else:                                                self.discarded_wrong_start += int(N_wrong_start)
        elif replace:                                            self.discarded_wrong_start = 'unknown'
        if N_no_cassette is not None:
            if replace or self.discarded_no_cassette=='unknown': self.discarded_no_cassette = int(N_no_cassette)
            else:                                                self.discarded_no_cassette += int(N_no_cassette)
        elif replace:                                            self.discarded_no_cassette = 'unknown'
        # special case for when we don't know the specific discarded categories, but we know total discarded is 0, 
        #  so the specific categories must be 0 too:
        if self.discarded_read_count == 0:   self.discarded_wrong_start, self.discarded_no_cassette = 0, 0

    def add_nonaligned_reads(self, N_all_non_aligned, N_unaligned, N_multiple_aligned, replace=False):
        """ Add not-None arg values to non_aligned_read_count, unaligned and multiple_aligned (or replace them).
        
        If the original values are 'unknown', or replace is True, replace instead of adding.
        If any of the args is None, don't modify the original value, unles replace is True, then set to 'unknown'.
        """
        if N_all_non_aligned is not None:
            if replace or self.unaligned=='unknown': self.non_aligned_read_count = int(N_all_non_aligned)
            else:                                    self.non_aligned_read_count += int(N_all_non_aligned)
        elif replace:                                self.non_aligned_read_count = 'unknown'
        if N_unaligned is not None:
            if replace or self.unaligned=='unknown': self.unaligned = int(N_unaligned)
            else:                                    self.unaligned += int(N_unaligned)
        elif replace:                                self.unaligned = 'unknown'
        if N_multiple_aligned is not None:
            if replace or self.multiple_aligned=='unknown': self.multiple_aligned = int(N_multiple_aligned)
            else:                                           self.multiple_aligned += int(N_multiple_aligned)
        elif replace:                                       self.multiple_aligned = 'unknown'
        # special case for when we don't know the specific unaligned categories, but we know total non-aligned is 0, 
        #  so the specific categories must be 0 too:
        if self.non_aligned_read_count==0:  self.unaligned, self.multiple_aligned = 0, 0

    @property
    def aligned_read_count(self):
        return sum([m.read_info(self.dataset_name).total_read_count for m in self.dataset])
    @property
    def perfect_read_count(self):
        return sum([m.read_info(self.dataset_name).perfect_read_count for m in self.dataset])
    @property
    def aligned_incl_removed(self):
        return self.aligned_read_count + sum(self.ignored_region_read_counts.values())

    @property
    def processed_read_count(self):
        """ Total processed readcount (integer): aligned + unaligned (ignored if unknown) + removed due to region. """
        known_values = self.aligned_read_count + sum(self.ignored_region_read_counts.values())
        try:                return known_values + self.non_aligned_read_count
        except TypeError:   return known_values
    @property
    def processed_read_count_str(self):
        """ Total processed readcount (string): aligned + unaligned ('unknown' if unknown) + removed due to region. """
        known_values = self.aligned_read_count + sum(self.ignored_region_read_counts.values())
        try:                return str(known_values + self.non_aligned_read_count)
        except TypeError:   return str(known_values) + '+unknown'

    @property
    def full_read_count(self):
        """ Full read count (integer): processed+discarded, or just processed if discarded is unknown. """
        try:                return self.processed_read_count + self.discarded_read_count
        except TypeError:   return self.processed_read_count
    @property
    def full_read_count_str(self):
        """ Full read count as a string: processed+discarded, or processed+'unknown' if discarded is unknown. """
        if 'unknown' not in self.processed_read_count_str:
            try:                return "%s"%(self.processed_read_count + self.discarded_read_count)
            except TypeError:   return "%s+unknown"%self.processed_read_count
        else:
            try:                return "%s+unknown"%(self.processed_read_count + self.discarded_read_count)
            except TypeError:   return "%s+unknown+unknown"%self.processed_read_count

    @property
    def strand_read_counts(self):
        strand_dict = {'+': 0, '-': 0}
        for m in self.dataset:
            if m.position.strand == 'both':
                for strand,count in m.original_strand_readcounts.items():
                    strand_dict[strand] += count
            else:
                strand_dict[m.position.strand] += m.read_info(self.dataset_name).total_read_count
        return strand_dict

    def reads_in_chromosome(self, chromosome):
        """ Return total number of reads in given chromosome."""
        return sum(m.read_info(self.dataset_name).read_info(self.dataset_name).total_read_count 
                   for m in self.dataset if m.position.chromosome==chromosome)

    @property
    def all_chromosomes(self):
        return set(mutant.position.chromosome for mutant in self.dataset)
    @property
    def cassette_chromosomes(self):
        return set(chrom for chrom in self.all_chromosomes if is_cassette_chromosome(chrom))
    @property
    def other_chromosomes(self):
        return set(chrom for chrom in self.all_chromosomes if is_other_chromosome(chrom))
    @property
    def non_genome_chromosomes(self):
        return self.cassette_chromosomes | self.other_chromosomes

    @property
    def mutants_in_genes(self):
        return len([1 for m in self.dataset if m.gene not in SPECIAL_GENE_CODES.all_codes])
    @property
    def mutants_not_in_genes(self):
        return len([1 for m in self.dataset if m.gene==SPECIAL_GENE_CODES.not_found])
    @property
    def mutants_undetermined(self):
        return len([1 for m in self.dataset if m.gene==SPECIAL_GENE_CODES.chromosome_not_in_reference])
    @property
    def mutant_counts_by_orientation(self):
        orientation_dict = defaultdict(int)
        for m in self.dataset:
            orientation_dict[m.orientation] += 1
        if '?' in orientation_dict:     del orientation_dict['?']
        if '-' in orientation_dict:     del orientation_dict['-']
        return orientation_dict
    @property
    def mutant_counts_by_feature(self):
        feature_dict = defaultdict(int)
        for m in self.dataset:
            feature_dict[m.gene_feature] += 1
        if '?' in feature_dict:  del feature_dict['?']
        if '-' in feature_dict:  del feature_dict['-']
        return feature_dict

    def mutants_in_chromosome(self, chromosome):
        """ Return total number of mutants in given chromosome."""
        return sum(1 for m in self.dataset 
                   if m.read_info(self.dataset_name).total_read_count and m.position.chromosome==chromosome)

    def merged_gene_feature_counts(self, merge_boundary_features=True, merge_confusing_features=False):
        """ Return (gene_feature,count) list, biologically sorted, optionally with all "boundary" features counted as one.

        The source gene feature counts are based on the self.mutant_counts_by_feature dict.
        If merge_confusing_features==True, any locations containing '??' will be listed as '??'.
        If merge_boundary_features==True, any locations containing '/' and no '??' will be listed as 'boundary'.
        The custom sort order (based on what seems sensible biologically) is: CDS, intron, UTR, other, boundary.
        """
        merged_feature_count_dict = defaultdict(int)
        for feature, count in self.mutant_counts_by_feature.items():
            # note that anything containing '??' AND '/' never gets merged as boundary
            if '??' in feature:
                if merge_confusing_features:                  merged_feature_count_dict['??'] += count
                else:                                         merged_feature_count_dict[feature] += count
            elif '/' in feature and merge_boundary_features:  merged_feature_count_dict['boundary'] += count
            else:                                             merged_feature_count_dict[feature] += count
        return merged_feature_count_dict

    @property
    def most_common_mutants(self):
        """ Return list of mutants with the most total reads (in dataset if multi-dataset)."""
        highest_readcount = max([mutant.read_info(self.dataset_name).total_read_count for mutant in self.dataset])
        highest_readcount_mutants = [mutant for mutant in self.dataset 
                                     if mutant.read_info(self.dataset_name).total_read_count==highest_readcount]
        return highest_readcount_mutants


class Insertional_mutant_pool_dataset():
    """ A dataset of insertional mutants - contains an iterable of Insertional_mutant objects, and a lot of extra data.

    May be a multi-dataset (i.e. a single dataset that actually contains readcounts for multiple samples) - in that case
     the Insertional_mutant objects will be multi-dataset mutants.

    WILL ADD MORE INFORMATION HERE ONCE I STOP CHANGING THIS SO OFTEN.
    
    Attributes - THIS MAY BE OUT OF DATE
     - cassette_end - specifies which end of the insertion cassette the reads are on, and 
     - reads_are_reverse - True if the reads are in reverse orientation to the cassette, False otherwise
     - discarded_read_count - number of reads discarded in preprocessing before alignment (not counted in processed_read_count)
     - ignored_region_read_counts - region_name:read_count dictionary (not counted in processed_read_count) ____ REALLY??
     - processed_read_count, unaligned_read_count - various read counts, obvious
     - strand_read_counts - name:count dictionaries to keep track of reads per strand
     - mutants_in_genes, mutants_not_in_genes, mutants_undetermined - counts of mutants in genes, not in genes, unknown
     - mutant_counts_by_orientation, mutant_count_by_feature - name:count dictionaries for mutant gene location properties

    For methods see method docstrings.
    """
    # TODO update docstring to contain up-to-date sensible info on everything!
    # - make sure to mention that all mutant positions are immutable by default, and how to deal with changing them
    # - explain about normal datasets and multi-datasets
    # - ____

    # Implement new functionality for datasets:
    # - MAYBE-TODO splitting joint datasets into separate ones?  Do we need that? PROBABLY NOT.
    # - MAYBE-TODO adding and subtracting datasets (by readcount) - do we currently need that?
    # - TODO-NEXT calculating growth rates!  Right now it's in mutant_growth_rates.py separately, but that may be temporary

    # MAYBE-TODO should I make multi-datasets a subclass of normal ones instead of having the same class implement both?

    def __init__(self, cassette_end='?', reads_are_reverse='?', multi_dataset=False, infile=None):
        """ Initializes empty dataset; saves properties as provided; optionally reads data from mutant infile. """
        # _mutants_by_position is the main data structure here, but it's also private:
        #      see "METHODS FOR EMULATING A CONTAINER TYPE" section below for proper ways to interact with mutants.
        #   the single mutants should be single-dataset by default, or multi-dataset if the containing dataset is.
        if multi_dataset:   blank_mutant_function = blank_multi_mutant_with_pos
        else:               blank_mutant_function = blank_single_mutant_with_pos
        self._mutants_by_position = keybased_defaultdict(blank_mutant_function)
        # various dataset summary data - single for a single dataset, a dictionary for a multi-dataset object.
        self.multi_dataset = multi_dataset
        if not multi_dataset:   self.summary = Dataset_summary_data(self, cassette_end, reads_are_reverse, None)
        else:                   self.summary = {}
        # data that's NOT related to a particular dataset
        # gene/annotation-related information - LATER-TODO should this even be here, or somewhere else?
        self.gene_annotation_header = []
        self.total_genes_in_genome = 0
        # optionally read mutant data from infile
        if infile is not None:
            self.read_data_from_file(infile)

    ######### METHODS FOR EMULATING A CONTAINER TYPE (a dataset is essentially a container of mutants)
    # everything is currently based on the private self._mutants_by_position dictionary, but this may change.
    # compared to just using a public self.mutants_by_position dictionary, this approach is better because:
    #   - it's more flexible
    #   - it doesn't allow having self.mutants_by_position[pos1] = mutant(pos2), with mismatched positions
    # Note that all of this will work just fine for multi-datasets without changes, since multi-datasets simply contain
    #  the same dictionary but with multi-dataset mutants instead of single-dataset ones.

    def __len__(self):      return len(self._mutants_by_position)
    def __iter__(self):     return self._mutants_by_position.itervalues()
    # for multi-datasets, for iterating over only mutants with 1+ reads in a particular dataset, see mutants_in_dataset()

    @property
    def size(self):         return len(self)

    # instead of __setitem__
    def add_mutant(self, mutant, overwrite=False):
        """ Add mutant to dataset. """
        if mutant.position in self._mutants_by_position.keys() and not overwrite:
            raise MutantError("Can't add mutant that would overwrite previous mutant at same position! "
                              +"Pass overwrite=True argument if you want to overwrite.")
        self._mutants_by_position[mutant.position] = mutant

    # instead of __delitem__
    def remove_mutant(self, *args, **kwargs):
        """ Remove mutant (by position) - can take a mutant, a position, or arguments to create a position. """
        if len(args)==1 and not kwargs and isinstance(args[0], Insertional_mutant):
            position = args[0].position
        else:
            position = self._make_position_object_if_needed(*args, **kwargs)
        del self._mutants_by_position[position]

    # instead of __getitem__
    def get_mutant(self, *args, **kwargs):
        """ Return the mutant with given position (given an Insertion_position instance, or arguments to create one).
        If mutant doesn't exist, create a new one with no reads/sequences. """
        return self._mutants_by_position[self._make_position_object_if_needed(*args,**kwargs)]

    # TODO how should mutant lookup by position deal with both-strand mutants?  You can have two mutants in the same position on opposite strands, but you CAN'T have one on both strands and one on + or -...  There should be safeguards to fold the + or - into the both-strand one if it's searched for, and things! THINK ABOUT THAT.

    # MAYBE-TODO implement get_Nth_mutant_by_position and get_Nth_mutant_by_readcount or some such?

    def __contains__(self, *args, **kwargs):
        """ Check if dataset contains mutant with given position (Insertion_position instance, or arguments to create one).

        If there's one (non-keyword) argument and it's an Insertion_position instance, use that for the check, 
         otherwise create a new Insertion_position instance with args/kwargs and check for that.

        You can only check "position in dataset" at present, not "mutant in dataset", since the latter would probably
         also just check by position, so the syntax would be misleading.
        """
        return self._make_position_object_if_needed(*args,**kwargs) in self._mutants_by_position
    
    @staticmethod
    def _make_position_object_if_needed(*args, **kwargs):
        """ Given either an Insertion_position instance or arguments to create one, return Insertion_position instance."""
        if len(args)==1 and not kwargs and isinstance(args[0], Insertion_position):
            return args[0]
        else:
            if 'immutable' not in kwargs:
                kwargs['immutable'] = True
            return Insertion_position(*args,**kwargs)


    ######### READING BASIC DATA INTO DATASET

    def add_alignment_reader_to_data(self, HTSeq_alignment_reader, uncollapse_read_counts=False, 
                                     treat_unknown_as_match=False, count_cassette=True, ignore_cassette=False):
        """ Adds all alignments from the reader to the mutant data; currently based only on position, but that may change. 

        Input must be a list/generator/etc of HTSeq.Alignment objects (usually an HTSeq.SAM_Reader).
        Set uncollapse_read_counts to True if the original deepseq data was collapsed to unique sequences using
         fastx_uncollapser before alignment, to get the correct original read counts.
        Treat_unknown_as_match governs whether alignments with no detailed information are treated as perfect or not.
        If count_cassette=True, add total cassette read count to header; 
         if ignore_cassette=True, ignore cassette reads in the data and list them as removed in the header.
        """
        if self.multi_dataset:  raise MutantError("add_alignment_reader_to_data not implemented for multi-datasets!")

        summ = self.summary
        if summ.cassette_end == '?':
            raise MutantError("Cannot add data from an alignment reader if cassette_end isn't specified! Please set the "
          +"summary.cassette_end attribute of this Insertional_mutant_pool_dataset instance to one of %s first."%SEQ_ENDS)
        if summ.reads_are_reverse == '?':
            raise MutantError("Cannot add data from an alignment reader if reads_are_reverse isn't set! Please set the "
                  +"reads_are_reverse attribute of this Insertional_mutant_pool_dataset instance to True/False first.")

        # Go over all the reads in the HTSeq_alignment_reader, add them to dataset
        for aln in HTSeq_alignment_reader:
            if uncollapse_read_counts:      read_count = get_seq_count_from_collapsed_header(aln.read.name)
            else:                           read_count = 1
            # if read is unaligned, add to unaligned count and skip to the next read
            if (not aln.aligned) or (aln.iv is None):
                summ.non_aligned_read_count += read_count
                continue
            # get the cassette insertion position (as an Insertion_position object) - USE IMMUTABLE POSITIONS BY DEFAULT
            position = get_insertion_pos_from_HTSeq_read_pos(aln.iv, summ.cassette_end, summ.reads_are_reverse, 
                                                             immutable_position=True)
            # if read is aligned to cassette and should be ignored, add to the right count and skip to the next read
            if ignore_cassette and is_cassette_chromosome(position.chromosome):
                summ.ignored_region_read_counts[position.chromosome] += read_count
                continue
            # MAYBE-TODO should summ.strand_read_counts just be generated from the data on the fly instead?
            # grab the right mutant based on the position, and add the reads to it; 
            curr_mutant = self.get_mutant(position)
            curr_mutant.add_read(aln, read_count, treat_unknown_as_match=treat_unknown_as_match)
        # special case for when we don't know the specific unaligned categories, but we know total non-aligned is 0, 
        #  so the specific categories must be 0 too:
        if summ.non_aligned_read_count==0:  summ.unaligned, summ.multiple_aligned = 0, 0

    def read_data_from_file(self, infile, assume_new_sequences=False):  # DEPRECATED!
        """ Read data from a file made by self.print_data, add mutants to dataset. Ignores some things. DEPRECATED. 

        Populates most of the dataset total read/mutant count values correctly, but ignores unaligned and discarded reads, 
         anything involving special regions, mutant-merging information, and possibly other things.
        Cannot deal with gene annotation fields.
        Ignores single sequence/count fields; the total number of sequence variants is unreliable if you add to preexisting
        data (if the data originally listed 1 unique sequence, and new data adds another 2 sequences, is the total 2 or 3 
         unique sequences?  If assume_new_sequences is True, the total is old+new; if it's False, it's max(old,new)). 

        DEPRECATED: All datasets should now have pickled versions for easy reading into python - use those instead.  
         This method will be kept around for reading old datasets without pickled versions, but will NOT BE UPDATED.
        """
        if self.multi_dataset:  raise MutantError("read_data_from_file not implemented for multi-datasets!")
        for line in open(infile):
            # try to extract info from the comment lines that can't be reproduced from the data, ignore the rest
            # NOTE: won't be adding any missing information - use pickled files instead, THIS IS DEPRECATED.
            if line.startswith('#'):                                        
                line = line.strip('#').strip()
                if line.startswith("Reads discarded in preprocessing"):
                    data = line.split('\t')[-1].split(' ')[0]
                    try:                data = int(data)
                    except ValueError:  pass
                    self.summary.discarded_read_count = data
                if line.startswith("Unaligned reads") or line.startswith("Reads without a unique alignment"):
                    data = line.split('\t')[-1].split(' ')[0]
                    try:                data = int(data)
                    except ValueError:  pass
                    self.summary.non_aligned_read_count = data
                    # special case for when we don't know the specific unaligned categories, but we know 
                    #  total non-aligned is 0, so the specific categories must be 0 too:
                    if self.summary.non_aligned_read_count==0:  self.summary.unaligned, self.summary.multiple_aligned = 0,0
                if line.startswith("(read location with respect to cassette: which end, which direction)"):
                    data = line.split('\t')[-1].strip('()').split(', ')
                    self.summary.cassette_end = data[0]
                    self.summary.reads_are_reverse = {'reverse':True, 'forward':False, '?':'?'}[data[1]]
                if line.startswith("(total genes in genome annotation data)"):
                    data = line.split('\t')[-1].strip('()')
                    try:                data = int(data)
                    except ValueError:  pass
                    self.total_genes_in_genome = data
                continue
            # ignore special-comment and header lines, parse other tab-separated lines into values
            if line.startswith('<REGEX>#') or line.startswith('<IGNORE>'):  continue       
            if line.startswith('chromosome\tstrand\tmin_position\t'):       continue
            # parse non-comment tab-separated lines into values
            fields = line.split('\t')
            chromosome = fields[0]
            strand = fields[1]
            min_pos = int(fields[2])
            full_pos = fields[3]
            gene, orientation, gene_feature = fields[4:7]
            total_reads,perfect_reads,sequence_variants = [int(x) for x in fields[7:10]]
            # generate new mutant if necessary; add counts and gene info to mutant (USE IMMUTABLE POSITIONS BY DEFAULT)
            position = Insertion_position(chromosome, strand, full_position=full_pos, immutable=True)
            curr_mutant = self.get_mutant(position)
            curr_mutant.add_counts(total_reads,perfect_reads,sequence_variants,assume_new_sequences)
            curr_mutant.update_gene_info(gene, orientation, gene_feature)
            # get however many specific sequences/counts are listed (this is variable)
            sequence_fields = fields[10::2]
            count_fields = fields[11::2]
            for seq, count in zip(sequence_fields, count_fields):
                if int(count)>0:
                    assert seq!=''
                    curr_mutant.sequences_and_counts[seq] += int(count)
            # add to dataset total read/mutant counts
            summ = self.summary
    

    ######### SUMMARY INFORMATION

    @property
    def mutants_by_gene(self):
        """ Return gene_name:mutant_list dict based on full list of dataset mutants; ignore SPECIAL_GENE_CODES) """
        mutants_by_gene = defaultdict(list)
        for mutant in self:
            if mutant.gene not in SPECIAL_GENE_CODES.all_codes:
                mutants_by_gene[mutant.gene].append(mutant)
        return mutants_by_gene
        # LATER-TODO add unit-test

    def get_gene_dict_by_mutant_number(self, dataset=None):
        """ Return mutant_count:gene_ID_set dict (genes with 1/2/etc mutants - in a particular dataset, if multi-dataset). 
        
        If the object is multi-dataset, dataset name must be provided, otherwise it cannot.
        """
        gene_by_mutantN = defaultdict(set)
        for (gene,mutants) in self.mutants_by_gene.iteritems():
            gene_by_mutantN[len([m for m in mutants if m.read_info(dataset).total_read_count])].add(gene)
        # check that the numbers add up to the total number of genes
        all_genes = set([mutant.gene for mutant in self]) - set(SPECIAL_GENE_CODES.all_codes)
        assert sum([len(geneset) for geneset in gene_by_mutantN.itervalues()]) == len(all_genes)
        return gene_by_mutantN
        #LATER-TODO add unit-test


    ######### MULTI-DATASET METHODS

    def populate_multi_dataset(self, source_dataset_dict, overwrite=False, check_gene_data=True):
        """ Given a dataset_name:single_dataset_object dictionary, populate current multi-dataset with the data. 

        If dataset already has a dataset with the same name, raise MutantError unless overwrite=True is passed.
        When merging mutants, don't double-check that their positions/genes match, unless check_gene_data=True is passed.
        Separate copies of mutant data are made, but the summary data object is added by reference.
        """
        if not self.multi_dataset:  raise MutantError("populate_multi_dataset can only be run for multi-datasets!")

        for dataset_name,dataset_object in source_dataset_dict.iteritems():

            if any([s in dataset_name for s in ' \t\n']):
                raise MutantError("Dataset name '%s' contains spaces/tabs/newlines - not allowed!"%dataset_name)
            if dataset_name in self.summary and not overwrite:
                raise MutantError("Dataset %s is already present in this multi-dataset! "%dataset_name
                                  +"Pass overwrite=True argument if you want to overwrite the previous data.")

            # copy source dataset summaries to own summary dict:
            #  (just doing a shallow-copy to avoid copying the whole dataset... copying ignored_region_read_counts
            summary_copy = copy.copy(dataset_object.summary)
            summary_copy.ignored_region_read_counts = defaultdict(int, dataset_object.summary.ignored_region_read_counts)
            summary_copy.dataset = self
            summary_copy.dataset_name = dataset_name
            self.summary[dataset_name] = summary_copy

            # Merge source dataset mutant data into own multi-dataset mutants
            #  (get_mutant will create new empty mutant if one doesn't exist)
            for mutant in dataset_object:
                curr_mutant = self.get_mutant(mutant.position)
                curr_mutant.add_other_mutant_as_dataset(mutant, dataset_name, overwrite=overwrite, 
                                                        check_constant_data=check_gene_data)
        # MAYBE-TODO add option to only include mutants with non-zero reads (total or perfect) in all datasets?  Or should that only happen during printing?  Or do we even care?  If I ever want to do that, there was code for it in the old version of mutant_join_datasets.py (before 2012-04-26)
        # This has no unit-tests, but has a run-test in mutant_join_datasets.py

    def mutants_in_dataset(self, dataset_name=None):
        """ List of all mutants with non-zero reads in dataset_name (or all mutants if dataset_name=None)."""
        return [mutant for mutant in self if dataset_name is None or mutant.by_dataset[dataset_name].total_read_count>0]

    def _check_dataset_consistency(self):
        """ Raise MutantError if self.summary, mutants, and sef.dataset_order don't all have the same set of datasets! """
        if not self.multi_dataset:  
            raise MutantError("_check_dataset_consistency only makes sense for multi-datasets!")
        def _check_sets_raise_error(set1, set2, set1_name, set2_name):
            if not set1==set2:
                raise MutantError("Multi-dataset mutant pool has different %s and %s dataset sets! %s, %s"%(set1_name, 
                                                                                              set2_name, set1, set2))
        datasets_from_summary = set(self.summary.keys())
        datasets_from_mutants = set.union(*[set(m.by_dataset.keys()) for m in self])
        _check_sets_raise_error(datasets_from_summary, datasets_from_mutants, "from summary", "from mutants")
        if self._dataset_order is not None:     
            datasets_from_order = set(self._dataset_order)
            _check_sets_raise_error(datasets_from_order, datasets_from_summary, "from dataset_order", "from summary")

    @property
    def dataset_order(self):
        """ A specific order of datasets, for printing - can be set directly, defaults to alphabetical sort. """
        self._check_dataset_consistency()
        try:                    return self._dataset_order
        except AttributeError:  return sorted(self.summary.keys())

    @dataset_order.setter
    def dataset_order(self, order):
        self._dataset_order = order
        self._check_dataset_consistency()


    ######### PROCESSING/MODIFYING DATASET

    def remove_mutants_based_on_other_dataset(self, other_dataset, readcount_min=1, perfect_reads=False):
        """ Remove any mutants with at least readcount_min reads in other_dataset (or perfect reads, if perfect_reads=True)

        This is based on EXACT position equality: a ?-101 mutant won't be removed due to a 100-? or 100-101 or 100-102 one.
        """
        # TODO do I want this to be based on non-exact position equality instead?
        if perfect_reads:   get_readcount = lambda m: m.perfect_read_count
        else:               get_readcount = lambda m: m.total_read_count
        # go over all mutants in self; need to convert the iterator to a list to make a separate copy, 
        #  otherwise we'd be modifying the iterator while iterating through it, which isn't allowed.
        for mutant in list(iter(self)):
            if get_readcount(other_dataset.get_mutant(mutant.position)) >= readcount_min:
                self.remove_mutant(mutant.position)
        # TODO really I shouldn't be removing mutants outright, just noting them as removed or something...  In that case should they or should they not show up in "for m in self"?  Probably not - they should have a separate dictionary?
        # TODO should I keep track of removed reads, and print in summary?  PROBABLY.

    def find_genes_for_mutants(self, genefile, detailed_features=False, N_run_groups=3, verbosity_level=1):
        """ To each mutant in the dataset, add the gene it's in (look up gene positions for each mutant using genefile).

        If detailed_features is True, also look up whether the mutant is in an exon/intron/UTR (NOT IMPLEMENTED); 
        Read the file in N_run_groups passes to avoid using up too much memory/CPU.
        """ 
        if self.multi_dataset:  raise MutantError("find_genes_for_mutants not implemented for multi-datasets!")
        # MAYBE-TODO implement for multi-datasets?  The actual gene-finding would be easy, since it'd just work on 
        #  multi-dataset mutants instead of single-dataset ones; adding stuff to summary would be harder.

        # group all the mutants by chromosome, so that I can go over each chromosome in genefile separately
        #   instead of reading in all the data at once (which uses a lot of memory)
        mutants_by_chromosome = defaultdict(set)
        for mutant in self:
            mutants_by_chromosome[mutant.position.chromosome].add(mutant)

        # First get the list of all chromosomes in the file, WITHOUT reading it all into memory
        with open(genefile) as GENEFILE:
            GFF_limit_data = GFF.GFFExaminer().available_limits(GENEFILE)
            chromosomes_and_counts = dict([(c,n) for ((c,),n) in GFF_limit_data['gff_id'].iteritems()])
            all_reference_chromosomes = set(chromosomes_and_counts.keys())

        # Now lump the chromosomes into N_run_groups sets with the feature counts balanced between sets, 
        #  to avoid using too much memory (by reading the whole file at once), 
        #   or using too much time (by reading the whole file for each chromosome/scaffold)
        chromosome_sets = split_into_N_sets_by_counts(chromosomes_and_counts, N_run_groups)

        ### go over all mutants on each chromosome, figure out which gene they're in (if any), keep track of totals
        # keep track of all the mutant and reference chromosomes to catch chromosomes that are absent in reference
        summ = self.summary
        self.total_genes_in_genome = 0
        for chromosome_set in chromosome_sets:
            genefile_parsing_limits = {'gff_id': list(chromosome_set)}
            if not detailed_features: 
                genefile_parsing_limits['gff_type'] = ['gene']
            with open(genefile) as GENEFILE:
                for chromosome_record in GFF.parse(GENEFILE, limit_info=genefile_parsing_limits):
                    if verbosity_level>1:   print "    parsing %s for mutant gene locations..."%chromosome_record.id
                    self.total_genes_in_genome += len(chromosome_record.features)
                    for mutant in mutants_by_chromosome[chromosome_record.id]:
                        gene_ID, orientation, feature = find_gene_by_pos(mutant.position, chromosome_record, 
                                                                         detailed_features, quiet=(verbosity_level==0))
                        mutant.gene, mutant.orientation, mutant.gene_feature = gene_ID, orientation, feature
                    if verbosity_level>1:   print "    ...found total %s genes."%(len(chromosome_record.features))
        if verbosity_level>1:   print "    found total %s genes in full genome."%(self.total_genes_in_genome)

        # for mutants in chromosomes that weren't listed in the genefile, use special values
        for chromosome in set(mutants_by_chromosome.keys())-set(all_reference_chromosomes):
            if not is_cassette_chromosome(chromosome):
                print 'Warning: chromosome "%s" not found in genefile data!'%(chromosome)
            for mutant in mutants_by_chromosome[chromosome]:
                mutant.gene,mutant.orientation,mutant.gene_feature = SPECIAL_GENE_CODES.chromosome_not_in_reference,'-','-'

    def add_gene_annotation(self, annotation_file, if_standard_Cre_file=False, custom_header=None, print_info=False):
        """ Add gene annotation to each mutant, based on annotation_file. See parse_gene_annotation_file doc for detail."""
        from parse_annotation_file import parse_gene_annotation_file
        gene_annotation_dict, gene_annotation_header = parse_gene_annotation_file(annotation_file, 
                     standard_Cre_file=if_standard_Cre_file, header_fields=custom_header, verbosity_level=int(print_info))
        # store the annotation header in self.summary, for printing
        if gene_annotation_header:  self.gene_annotation_header = gene_annotation_header
        else:                       self.gene_annotation_header = 'GENE_ANNOTATION_DATA'
        # add the annotation info to each mutant (or nothing, if gene has no annotation)
        # MAYBE-TODO should I even store gene annotation in each mutant, or just keep a separate per-gene dictionary?
        for mutant in self:
            try:                mutant.gene_annotation = gene_annotation_dict[mutant.gene]
            except KeyError:    mutant.gene_annotation = []
        # LATER-TODO add this to the gene-info run-test case!


    ######### MUTANT-MERGING AND ADJACENT-COUNTING

    def _possibly_adjacent_positions(self, max_distance, same_strand_only=False, include_cassette_chromosomes=True, 
                                     include_other_chromosomes=False):
        """ Generates all mutant position pairs that might be within max_distance of each other (same strand, or not).

        May or may not be well optimized (the most naive way would just be combinations(all_positions,2)).

        If same_strand_only is True, both-strand positions are counted as both +strand and -strand.
        """
        for chromosome in self.summary.all_chromosomes:
            if include_cassette_chromosomes==False and is_cassette_chromosome(chromosome): continue
            if include_other_chromosomes==False and is_other_chromosome(chromosome):       continue

            all_positions = sorted([mutant.position for mutant in self if mutant.position.chromosome==chromosome])

            # MAYBE-TODO if same_strand_only, could split up by strand, too!  (both-strand ones should go on both lists?)
            # MAYBE-TODO could optimize this further by only considering pairs of close-by positions, instead of all pairs:
            #  positions are sorted by min_position and then strand, so if you sorted all_positions, 
            #  two mutants on the same strand X bp away can be at most 2X positions away on the sorted position list, 
            #  so use a sliding window of size 2X+1.  #  (or of size X+1 if they're all on the same strand already)

            # go over each pos1,pos2 pair (only once)
            for (pos1, pos2) in combinations(all_positions,2):
                yield (pos1, pos2)

    def _merge_two_mutants(self, pos1, pos2):
        """ Merge two mutants (given by position), while making sure not to mess up dictionary key immutability. 
        
        Doing the merge may require the position to change, and changing a dictionary key is BAD, 
         so first remove both mutants from the dataset, make the position mutable, 
         merge the mutants (using the Insertional_mutant.merge_mutant method, which does all the real work), 
         make the position immutable again, and finally add the merged mutant back to the dataset.
        """
        mutant1 = self.get_mutant(pos1)
        mutant2 = self.get_mutant(pos2)
        self.remove_mutant(pos1)
        self.remove_mutant(pos2)
        mutant1.position.make_mutable_REMEMBER_CLEANUP_FIRST()
        mutant1.merge_mutant(mutant2, opposite_strand_tandem=(pos1.strand!=pos2.strand))
        mutant1.position.make_immutable()
        self.add_mutant(mutant1)

    def merge_adjacent_mutants(self, merge_max_distance=1, merge_count_ratio=1000, merge_cassette_chromosomes=False, 
                               merge_other_chromosomes=True, OUTPUT=sys.stdout):
        """ Merge adjacent mutants based on strand, distance, and count ratio; save counts.

        Merge mutants if they're distant by merge_max_distance or less and one has at least merge_count_ratio x fewer 
         reads than the other; merge their read counts, sequences etc, and remove the lower-read-count one from dataset.

        Cannot be run after merge_opposite_tandem_mutants - NOT IMPLEMENTED for both-stranded mutants.

        If merge_cassette_chromosomes or merge_other_chromosomes is False, don't do merging on cassette chromosomes 
         and other non-genome chromosomes (chloroplast, mitochondrial, etc) respectively.

        Write details to OUTPUT (open filehandle; pass sys.stdout if desired, or None for no writing);
         Save information on merged counts to 
        """
        if self.multi_dataset:  raise MutantError("merge_adjacent_mutants not implemented for multi-datasets!")
        if merge_count_ratio < 1:  raise MutantError("merge_count_ratio must be at least 1!")
        if self.summary.adjacent_merging_count_ratio is not None:
            raise MutantError("merge_adjacent_mutants has been run already - cannot run twice on same dataset!")
            # MAYBE-TODO or should that be allowed?  But then it could be done with different conditions, which would 
            #  make things confusing, and the new self.summary.adjacent_merging_count_ratio would overwrite the old...
        if self.summary.adjacent_max_distance not in (None, merge_max_distance):
            raise MutantError("Cannot change max adjacent distance once it's been set! (by merging or counting) "
                              "Currently %s, can't change to %s."%(self.summary.adjacent_max_distance, merge_max_distance))
            # MAYBE-TODO or allow this and just re-run count_adjacent_mutants afterward with new value?
        if not self.summary.merging_which_chromosomes in [(merge_cassette_chromosomes,merge_other_chromosomes), 
                                                          (None,None)]:
            raise MutantError("Cannot change which chromosomes are subject to mutant-merging! "
                             +"Currently %s for cassette, %s for non-nuclear."%self.summary.merging_which_chromosomes
                             +"Trying to change to %s and %s."%(merge_cassette_chromosomes, merge_other_chromosomes))
        if any([mutant.position.strand=='both' for mutant in self]):
            raise MutantError("merge_adjacent_mutants should be run BEFORE merge_opposite_tandem_mutants "
                             +"- NOT IMPLEMENTED for both-stranded mutants!")
            # LATER-TODO implement for both-stranded mutants too! In that case all strand-comparisons should be
            #  rewritten so 'both' and either + or - would compare as the same, or something... 
        if OUTPUT is None: OUTPUT = FAKE_OUTFILE
        OUTPUT.write("# Merging adjacent mutants: max distance %s, min count ratio %s\n"%(merge_max_distance, 
                                                                                          merge_count_ratio))
        OUTPUT.write(" (The merged mutant will have the position of whichever original mutant had more reads)\n")
        if merge_count_ratio == 1:  
            OUTPUT.write(" (Warning: minimum ratio is 1, so sometimes both mutants will have the same read number "
                         "- the earlier position is used in that case.)\n")
        OUTPUT.write(" (%sincluding cassette chromosomes)"%('' if merge_cassette_chromosomes else 'not ')
                    +" (%sincluding non-nuclear chromosomes)\n"%('' if merge_other_chromosomes else 'not '))
        merged_count = 0
        for pos1,pos2 in self._possibly_adjacent_positions(merge_max_distance, True, merge_cassette_chromosomes, 
                                                           merge_other_chromosomes):
            # if the two positions are on different chromosomes or strands or aren't adjacent, skip
            if pos1.chromosome != pos2.chromosome:                              continue
            if pos1.strand != pos2.strand:                                      continue
            if abs(pos2.min_position-pos1.min_position) > merge_max_distance:   continue
            # all other mutant pairs should be same-strand adjacent - up for merging depending on readcount ratio.

            # make sure these mutants still exist in the dataset (haven't been merged) (note: if they don't exist, 
            #   the dataset will silently create new mutants with readcounts 0, so check "x in self" explicitly)
            #  (this shouldn't be a problem with a situation like 1,1000,1 readcounts - the first and third will both
            #   be merged into the second, which won't disappear, so all is fine. However with a situation like
            #   1000,1,1000 readcounts, the middle mutant will be merged into only one of the sides, which is right.
            #   Which side it is will be arbitrary, but I think that's fine for now - MAYBE-TODO make it random, 
            #    or split the reads between the two, or use the earlier-position one? Randomness is hard to test...)
            if pos1 in self:
                mutant1 = self.get_mutant(pos1)
            else:
                OUTPUT.write("Warning: attempting to merge a mutant that no longer exists %s with %s!\n"%(pos1, pos2))
                continue
            if pos2 in self:
                mutant2 = self.get_mutant(pos2)
            else:
                OUTPUT.write("Warning: attempting to merge a mutant that no longer exists %s with %s!\n"%(pos2, pos1))
                continue

            # if the two mutants are adjacent and with sufficiently different readcounts, count, 
            #  merge the lower-readcount mutant into the higher-readcount one, and remove the lower one from the dataset
            mutant1_readcount, mutant2_readcount = mutant1.total_read_count, mutant2.total_read_count
            readcount_ratio = max(mutant1_readcount/mutant2_readcount, mutant2_readcount/mutant1_readcount)
            if not readcount_ratio>=merge_count_ratio:  continue
            merged_count += 1
            OUTPUT.write(" MERGING same-strand adjacent mutants: %s and %s, %s and %s reads.\n"%(pos1, pos2, 
                                                                                mutant1_readcount, mutant2_readcount))
            assert pos1.strand==pos2.strand, "Mutants for adjacent-merging must be on same strand!"
            # make sure that mutant1.total_read_count >= mutant.total_read_count2 - switch if they're otherwise
            if mutant1_readcount < mutant2_readcount:
                pos1, pos2 = pos2, pos1
            # if both have equal readcounts, use the earlier one (pos1)
            if mutant1_readcount == mutant2_readcount:
                OUTPUT.write("  (mutants have same readcount - using earlier position, %s)\n"%pos1)
            self._merge_two_mutants(pos1, pos2)

        # add overall counts to dataset summary
        OUTPUT.write("# Finished merging adjacent mutants: %s pairs merged\n"%merged_count) 
        self.summary.adjacent_max_distance = merge_max_distance
        self.summary.adjacent_merging_count_ratio = merge_count_ratio
        self.summary.merging_which_chromosomes = (merge_cassette_chromosomes, merge_other_chromosomes)
        self.summary.merged_adjacent_pairs = merged_count

    # TODO there was a bug here where mutants 1,2,3 all should be merged together and a warning gets printed!  See the warning in test_data/count-aln__merge-adjacent2-r3_merging-info.txt, which really shouldn't be there.  Research and hopefully fix that!  I'm not sure whether the problem is ONLY that is prints an unnecessary warning, or if there's also an actual issue with merging.
    # TODO was there some other bug when doing tandem-merging and adjacent-merging at once?  I think something weird came up in actual data analysis - see ../../1206_Ru-screen1_deepseq-data-early/notes.txt  "Mutants" section.
    # LATER-TODO should put a header line on merging-outfile giving the main outfile name; and also give the main outfile a header line saying that the mutant merging info is in options.mutant_merging_outfile.

    def merge_opposite_tandem_mutants(self, merge_cassette_chromosomes=False, merge_other_chromosomes=True, 
                                      OUTPUT=sys.stdout):
        """ Merge opposite-strand tandem mutants (in same position but opposite strands), setting strand to 'both'. 

        Merge mutants if they're in the same position but on opposite strands; merge their read counts, sequences etc, 
         remove both from the dataset, and add a new one to the dataset, with 'both' for strand and gene orientation.
         If the full positions of the two mutants were 100-? and ?-101, the new one is 100-101.

        If merge_cassette_chromosomes or merge_other_chromosomes is False, don't do merging on cassette chromosomes 
         and other non-genome chromosomes (chloroplast, mitochondrial, etc) respectively.

        Write details to OUTPUT (open filehandle; pass sys.stdout if desired, or None for no printing);
         save the number of merged cases to self.summary.
        """
        if self.multi_dataset:  raise MutantError("merge_opposite_tandem_mutants not implemented for multi-datasets!")
        if not self.summary.merging_which_chromosomes in [(merge_cassette_chromosomes,merge_other_chromosomes), 
                                                          (None,None)]:
            raise MutantError("Cannot change which chromosomes are subject to mutant-merging! "
                              "Currently %s for cassette, %s for non-nuclear."%self.summary.merging_which_chromosomes)
        if OUTPUT is None: OUTPUT = FAKE_OUTFILE
        OUTPUT.write("# Merging opposite-strand same-position mutants (presumably tail-to-tail tandems)\n"
                    +" (%sincluding cassette chromosomes)"%('' if merge_cassette_chromosomes else 'not ')
                    +" (%sincluding non-nuclear chromosomes)\n"%('' if merge_other_chromosomes else 'not '))

        merged_tandem_count = 0
        all_positions = sorted([mutant.position for mutant in self])
        # same-position opposite-strand mutants will always be adjacent on the sorted position list, so only look at those
        for pos1,pos2 in zip(all_positions,all_positions[1:]):
            # if the two mutants are on different chromosomes or different positions, go on to the next pair
            if pos1.chromosome != pos2.chromosome:        continue
            if pos2.min_position != pos1.min_position:    continue
            # if the two mutants are on chromosomes that shouldn't be merged, go on to the next pair
            if merge_cassette_chromosomes==False and is_cassette_chromosome(pos1.chromosome): continue
            if merge_other_chromosomes==False and is_other_chromosome(pos1.chromosome):       continue
            assert pos1.strand != pos2.strand, "Mutants with same chromosome/position must have different strands!"
            assert 'both' not in [pos1.strand, pos2.strand], "Two mutants in one position must be +/- strand, not both!"

            OUTPUT.write(" MERGING opposite-strand same-position tandem mutants: %s and %s.\n"%(pos1,pos2))
            merged_tandem_count += 1
            self._merge_two_mutants(pos1, pos2)

        # add overall counts to dataset summary
        OUTPUT.write("# Finished merging opposite-strand same-position mutants: %s pairs merged\n"%merged_tandem_count)
        self.summary.merging_which_chromosomes = (merge_cassette_chromosomes, merge_other_chromosomes)
        self.summary.merged_opposite_tandems = merged_tandem_count

    def count_adjacent_mutants(self, adjacent_max_distance=1, count_cassette_chromosomes=False, 
                               count_other_chromosomes=True, different_parameters=False, OUTPUT=sys.stdout):
        """ Count various categories of adjacent mutants; save data to self.summary. 
        
        Specific cases counted are:
            - adjacent mutants on the same strand
            - same-position mutants on opposite strands (100-? and ?-101 - both have min_position 100)
            - adjacent mutants on opposite strands, facing away from each other (100-? and ?-103)
            - adjacent mutants on opposite strands, facing toward each other (?-101 and 102-?)
            (the reason for differentiating between the "away" and "toward" cases is that 
             the "away" case could be a tail-to-tail tandem with a genomic deletion in the middle, 
             but the "toward" case really HAS TO be random unrelated mutants (right?), so it can be used as a reference.)

        If count_cassette_chromosomes or count_other_chromosomes is False, don't include in the count cassette chromosomes 
         and other non-nuclear chromosomes (chloroplast, mitochondrial, etc) respectively.
        If different_parameters is True, check that adjacent_max_distance, count_cassette_chromosomes, 
         and count_other_chromosomes match previously used merging settings.


        Also write details to OUTPUT (open filehandle; pass sys.stdout if desired, or None for no printing).
        """
        if self.multi_dataset:  raise MutantError("count_adjacent_mutants not implemented for multi-datasets!")
        if not different_parameters:
            if not self.summary.merging_which_chromosomes in [(count_cassette_chromosomes,count_other_chromosomes), 
                                                              (None,None)]:
                raise MutantError("If trying to change which chromosomes are included when doing adjacent-mutant-count "
                                  "compared to mutant-merging, use different_parameters argument."
                                  "(Merging: %s for cassette, %s for non-nuclear.)"%self.summary.merging_which_chromosomes)
            if not self.summary.adjacent_max_distance in [adjacent_max_distance, None]:
                raise MutantError("If trying to change the adjacent-mutant max distance when doing adjacent-mutant-count "
                                  "compared to mutant-merging, use different_parameters argument."
                                  "(Merging: max distance %s.)"%self.summary.adjacent_max_distance)
        if OUTPUT is None: OUTPUT = FAKE_OUTFILE
        OUTPUT.write("# Counting adjacent mutants: max distance %s\n"%adjacent_max_distance
                    +" (%sincluding cassette chromosomes)"%('' if count_cassette_chromosomes else 'not ')
                    +" (%sincluding non-nuclear chromosomes)\n"%('' if count_other_chromosomes else 'not '))
        adjacent_same_strand = 0
        same_position_opposite_strands = 0
        adjacent_opposite_strands_away = 0
        adjacent_opposite_strands_toward = 0
        # MAYBE-TODO may want to make all the adjacent_* counts dictionaries by distance?
        for pos1,pos2 in self._possibly_adjacent_positions(adjacent_max_distance, False, count_cassette_chromosomes, 
                                                           count_other_chromosomes):
            # if the two positions are on different chromosomes or aren't adjacent, skip
            if pos1.chromosome != pos2.chromosome:                                continue
            if abs(pos2.min_position-pos1.min_position) > adjacent_max_distance:  continue
            # same position, opposite strands
            if pos1.min_position==pos2.min_position:
                assert pos1.strand != pos2.strand, "Two mutants with same position and strand shouldn't happen!"
                assert 'both' not in (pos1.strand,pos2.strand), "A both-strand mutant can't be same-position to another!"
                same_position_opposite_strands += 1
                OUTPUT.write("  opposite-strand same-position tandem mutants: %s and %s.\n"%(pos1,pos2))
            # adjacent positions, same strand - remember that both-strand mutants are same-strand with everything!
            elif pos1.strand == pos2.strand or 'both' in (pos1.strand,pos2.strand): 
                assert pos1.min_position!=pos2.min_position, "Two mutants with same position and strand shouldn't happen!"
                adjacent_same_strand += 1
                mutant1_readcount = self.get_mutant(pos1).total_read_count
                mutant2_readcount = self.get_mutant(pos2).total_read_count
                OUTPUT.write("  same-strand adjacent mutants: %s and %s, %s and %s reads.\n"%(pos1,pos2, 
                                                                                     mutant1_readcount, mutant2_readcount))
            # adjacent positions, opposite strands
            else:
                assert pos1.min_position != pos2.min_position and pos1.strand != pos2.strand, "Not adjacent-opposite!"
                assert 'both' not in (pos1.strand,pos2.strand), "A both-strand mutant can't be opposite-strand to another!"
                #  ?-X is -strand, X-? is +strand.  
                # away  =  100-? and ?-103  =  +strand has a lower position than -strand;  toward  =  other way around
                first_pos = min([pos1, pos2])
                if first_pos.strand == '+':
                    OUTPUT.write('  adjacent opposite-strand "away-facing" mutants: %s and %s.\n'%(pos1,pos2))
                    adjacent_opposite_strands_away += 1
                else:
                    OUTPUT.write('  adjacent opposite-strand "toward-facing" mutants: %s and %s.\n'%(pos1,pos2))
                    adjacent_opposite_strands_toward += 1

        # add overall counts to dataset summary
        OUTPUT.write("# Finished counting adjacent mutants (max distance %s): "%adjacent_max_distance
                     +"%s adjacent same-strand pairs, "%adjacent_same_strand
                     +"%s same-position opposite-strand pairs, "%same_position_opposite_strands
                     +"%s adjacent opposite-strand \"away-facing\" pairs "%adjacent_opposite_strands_away
                     +"(may be tandems with a deletion), "
                     +"%s adjacent opposite-strand \"toward-facing\" pairs "%adjacent_opposite_strands_toward
                     +"(definitely two separate mutants).\n")
        self.summary.adjacent_same_strand = adjacent_same_strand
        self.summary.same_position_opposite = same_position_opposite_strands
        self.summary.adjacent_opposite_away = adjacent_opposite_strands_away
        self.summary.adjacent_opposite_toward = adjacent_opposite_strands_toward
        self.summary.adjacent_max_distance = adjacent_max_distance
        if not different_parameters:
            self.summary.merging_which_chromosomes = (count_cassette_chromosomes, count_other_chromosomes)


    ######### WRITING DATA TO FILES

    def _get_summary(self, dataset=None):
        """ Help function to unify single/multi-datasets: return summary if dataset=None, else summary[dataset]. """
        # TODO should this be in some other section?
        if dataset is None:     return self.summary
        else:                   return self.summary[dataset]

    def _most_common_mutants_info(self, dataset=None):
        """ Return a string containing details for most common mutant, or count of most common mutants if multiple. """
        summ = self._get_summary(dataset)
        most_common_mutants = summ.most_common_mutants
        m = most_common_mutants[0]
        # calculate the fraction of total reads per mutant, assuming each mutant has the same readcount
        assert len(set([m.read_info(dataset).total_read_count for m in most_common_mutants])) == 1
        readcount_info = value_and_percentages(m.read_info(dataset).total_read_count, [summ.aligned_read_count])
        if len(most_common_mutants) == 1:   return "%s (%s)"%(readcount_info, m.position)
        else:                               return "%s (%s mutants)"%(readcount_info, len(most_common_mutants))

    @staticmethod
    def _sort_feature_list(feature_list):
        """ Sort feature list: biologically (CDS, intron, UTR, other, boundary), then alphabetically in each category. """
        # using a dictionary to assign each category a sorting position, with "other" (default) between UTR and boundary.
        proper_feature_order = defaultdict(lambda: 3, 
                                           {'CDS':0, 'intron':1, 'five_prime_UTR':2, 'three_prime_UTR':2, 'boundary':4})
        feature_count_sorted_list = sorted(feature_list, key=lambda f: (proper_feature_order[f],f))
        return feature_count_sorted_list

    @staticmethod
    def _make_genelist_str(geneset, N_genes_to_print):
        """ Given a set of names, return a string containing at most N_genes of the names: '(A,B,...)', '(A,B)' or ''. """
        if N_genes_to_print>0:
            genelist_to_print = ', '.join(sorted(geneset)[:N_genes_to_print])
            if len(geneset)<=N_genes_to_print:  return '(%s)'%genelist_to_print
            else:                               return '(%s, ...)'%genelist_to_print
        else:                                   return ''

    def print_summary(self, OUTPUT=sys.stdout, N_genes_to_print=5, line_prefix='    ', header_prefix=' * ', 
                      merge_boundary_features=True, count_cassette=True, count_other=True):
        """ Print basic summary info about the dataset/s: read and mutant counts and categories, gene numbers, etc.

        Prints tab-separated table for multi-datasets.
        Prints to stdout by default, can also pass an open file object).
        """
        ### define a list of datasets+summaries that we'll be dealing with
        if not self.multi_dataset:  summaries_and_datasets = [(self.summary, None)]
        else:                  summaries_and_datasets = [(self.summary[dataset],dataset) for dataset in self.dataset_order]
        summaries, datasets = zip(*summaries_and_datasets)

        ### First set up a list of line descriptions and value-getter functions
        # Note on lambdas: wherever there's a for loop, I have to use lambda x=x (region,orientation,etc) 
        #  to bind x to the CURRENT value of x - otherwise whenever the lambda got called, 
        #   the current environment value of x would be used instead.
        descriptions_and_value_getters = DVG = []

        DVG.append((header_prefix+"Total reads in sample:", 
                    lambda summ,mutants,dataset: "%s"%(summ.full_read_count_str) ))
        def _fraction_or_unknown(value, totals):
            try:                return value_and_percentages(value, totals)
            except TypeError:   return "%s (unknown)"%value
        DVG.append((header_prefix+"Reads discarded in preprocessing (% of total):", 
                lambda summ,mutants,dataset: _fraction_or_unknown(summ.discarded_read_count, [summ.full_read_count])))
        DVG.append((line_prefix+"discarded due to wrong start (% of total):", 
                lambda summ,mutants,dataset: _fraction_or_unknown(summ.discarded_wrong_start, [summ.full_read_count])))
        DVG.append((line_prefix+"discarded due to no cassette (% of total):", 
                lambda summ,mutants,dataset: _fraction_or_unknown(summ.discarded_no_cassette, [summ.full_read_count])))

        DVG.append((header_prefix+"Reads without a unique alignment (% of total, % of post-preprocessing):", 
                    lambda summ,mutants,dataset: _fraction_or_unknown(summ.non_aligned_read_count, 
                                                               [summ.full_read_count, summ.processed_read_count]) ))
        DVG.append((line_prefix+"unaligned reads (% of total, % of post-preprocessing):", 
                    lambda summ,mutants,dataset: _fraction_or_unknown(summ.unaligned, 
                                                               [summ.full_read_count, summ.processed_read_count]) ))
        DVG.append((line_prefix+"multiply aligned reads (% of total, % of post-preprocessing):", 
                    lambda summ,mutants,dataset: _fraction_or_unknown(summ.multiple_aligned, 
                                                               [summ.full_read_count, summ.processed_read_count]) ))

        DVG.append((header_prefix+"Uniquely aligned reads (% of total, % of post-preprocessing):",
                    lambda summ,mutants,dataset: value_and_percentages(summ.aligned_incl_removed, 
                                                               [summ.full_read_count, summ.processed_read_count]) ))

        all_ignored_regions = set.union(*[set(summ.ignored_region_read_counts) for summ in summaries])
        for region in sorted(all_ignored_regions):
            DVG.append((line_prefix+"Removed reads aligned to %s (%% of total, %% of post-preprocessing):"%region, 
                        lambda summ,mutants,dataset,region=region: value_and_percentages(
                                                                summ.ignored_region_read_counts[region], 
                                                                [summ.full_read_count, summ.processed_read_count]) ))
        if all_ignored_regions:
            DVG.append((header_prefix+"Remaining aligned reads (% of total, % of post-preprocessing):", 
                        lambda summ,mutants,dataset: value_and_percentages(summ.aligned_read_count, 
                                                               [summ.full_read_count, summ.processed_read_count]) ))

        DVG.append((line_prefix+"Perfectly aligned reads, no mismatches (% of aligned):", 
                    lambda summ,mutants,dataset: value_and_percentages(summ.perfect_read_count, 
                                                                       [summ.aligned_read_count]) ))

        for strand in sorted(set.union(*[set(summ.strand_read_counts) for summ in summaries])):
            DVG.append((line_prefix+"Reads with cassette direction matching chromosome %s strand (%% of aligned):"%strand,
                        lambda summ,mutants,dataset,strand=strand: value_and_percentages(summ.strand_read_counts[strand], 
                                                                                         [summ.aligned_read_count]) ))

        special_chromosomes = []
        if count_cassette:
            special_chromosomes += sorted(summ.cassette_chromosomes)
        if count_other:
            special_chromosomes += sorted(summ.other_chromosomes)

        for chromosome in special_chromosomes:
            DVG.append((line_prefix+"Reads aligned to %s (%% of aligned):"%chromosome, 
                        lambda summ,mutants,dataset,chromosome=chromosome: value_and_percentages(
                                                                              summ.reads_in_chromosome(chromosome),
                                                                              [summ.aligned_read_count]) ))
        
        DVG.append((header_prefix+"Mutant merging/counts (deciding when different-position reads should be one mutant)", 
                    lambda summ,mutants,dataset: '' ))
        DVG.append((line_prefix+" (adjacent-merging/counting max distance):", 
                    lambda summ,mutants,dataset: "(%s)"%summ.adjacent_max_distance ))
        DVG.append((line_prefix+" (if we're including mutants in cassette and in non-nuclear chromosomes):", 
                    lambda summ,mutants,dataset: "(%s, %s)"%summ.merging_which_chromosomes )) 
        DVG.append((line_prefix+"merged same-strand adjacent mutant pairs and opposite-strand tandem pairs:", 
                    lambda summ,mutants,dataset: "%s, %s"%(summ.merged_adjacent_pairs, summ.merged_opposite_tandems) ))
        DVG.append((line_prefix+" (adjacent-merging min count ratio - None if no adjacent-merging):", 
                    lambda summ,mutants,dataset: "(%s)"%(summ.adjacent_merging_count_ratio) ))
        DVG.append((line_prefix+"remaining same-position opposite-strand pairs (if not merged as tandems):", 
                    lambda summ,mutants,dataset: "%s"%summ.same_position_opposite )) 
        DVG.append((line_prefix+'remaining adjacent opposite-strand "toward-facing" pairs (those are definitely real):', 
                    lambda summ,mutants,dataset: "%s"%summ.adjacent_opposite_toward )) 
        DVG.append((line_prefix+'remaining adjacent opposite-strand "away-facing" pairs (% of toward-facing):', 
                    lambda summ,mutants,dataset: value_and_percentages(summ.adjacent_opposite_away, 
                                                       [summ.adjacent_opposite_toward], percentage_format_str='%.0f') )) 
        DVG.append((line_prefix+'remaining adjacent same-strand unmerged pairs (% of 2*toward-facing):', 
                    lambda summ,mutants,dataset: value_and_percentages(summ.adjacent_same_strand, 
                                                       [2*summ.adjacent_opposite_toward], percentage_format_str='%.0f') )) 

        DVG.append((header_prefix+"Distinct mutants (read groups) by cassette insertion position:", 
                    lambda summ,mutants,dataset: "%s"%len(mutants) ))
        DVG.append((line_prefix+"(read location with respect to cassette: which end, which direction):", 
                    lambda summ,mutants,dataset: "(%s, %s)"%(summ.cassette_end, 
                                                {'?': '?', True: 'reverse', False: 'forward'}[summ.reads_are_reverse]) ))
        DVG.append((line_prefix+"(average and median reads per mutant):", 
                    lambda summ,mutants,dataset: "(%d, %d)"%(round((summ.aligned_read_count)/len(mutants)), 
                                             round(median([m.read_info(dataset).total_read_count for m in mutants]))) ))

        DVG.append((line_prefix+"Most common mutant(s): reads (% of aligned) (position or count):",
                    lambda summ,mutants,dataset: self._most_common_mutants_info(dataset) ))
        # MAYBE-TODO may also be a good idea to keep track of the most common SEQUENCE, not just mutant...

        for chromosome in special_chromosomes:
            DVG.append((line_prefix+"Mutant cassettes in %s (%% of total):"%chromosome, 
                        lambda summ,mutants,dataset,chromosome=chromosome: value_and_percentages(
                                                                            summ.mutants_in_chromosome(chromosome),
                                                                            [len(mutants)]) ))

        # print the gene annotation info, but only if there is any
        if any([summ.mutants_in_genes+summ.mutants_not_in_genes+summ.mutants_undetermined for summ in summaries]):
            DVG.append((line_prefix+"Mutant cassettes with unknown gene info (probably cassette-mapped) (% of total):", 
                        lambda summ,mutants,dataset: value_and_percentages(summ.mutants_undetermined, [len(mutants)]) ))
            DVG.append((line_prefix+"Mutant cassettes in intergenic spaces (% of total, % of known):", 
                        lambda summ,mutants,dataset: value_and_percentages(summ.mutants_not_in_genes, 
                                               [len(mutants), summ.mutants_not_in_genes+summ.mutants_in_genes]) ))
            DVG.append((header_prefix+"Mutant cassettes inside genes (% of total, % of known):", 
                        lambda summ,mutants,dataset: value_and_percentages(summ.mutants_in_genes, 
                                               [len(mutants), summ.mutants_not_in_genes+summ.mutants_in_genes]) ))
            for orientation in sorted(set.union(*[set(summ.mutant_counts_by_orientation) for summ in summaries]), 
                                      reverse=True):
                DVG.append((line_prefix+"Mutant cassettes in %s orientation to gene (%% of ones in genes):"%orientation, 
                            lambda summ,mutants,dataset,o=orientation: value_and_percentages(
                                                                                summ.mutant_counts_by_orientation[o], 
                                                                                [summ.mutants_in_genes]) ))
            # custom order for features to make it easier to read: CDS, intron, UTRs, everything else alphabetically after
            # MAYBE-TODO also give print_summary an option for merge_confusing_features arg to merged_gene_feature_counts?
            for feature in self._sort_feature_list(set.union(
                                *[set(summ.merged_gene_feature_counts(merge_boundary_features)) for summ in summaries])):
                DVG.append((line_prefix+"Mutant cassettes in gene feature %s (%% of ones in genes):"%feature, 
                            lambda summ,mutants,dataset,feature=feature: value_and_percentages(
                                                    summ.merged_gene_feature_counts(merge_boundary_features)[feature], 
                                                    [summ.mutants_in_genes]) ))

            _N_all_genes = lambda dataset: sum([len(genes) for N_mutants,genes 
                                               in self.get_gene_dict_by_mutant_number(dataset).items() if N_mutants>0])
            _N_genes_in_multiple_mutants = lambda dataset: sum([len(genes) for N_mutants,genes 
                                               in self.get_gene_dict_by_mutant_number(dataset).items() if N_mutants>1])
            DVG.append((header_prefix+"Genes containing a mutant (% of all genes):", 
                        lambda summ,mutants,dataset: value_and_percentages(_N_all_genes(dataset), 
                                                                           [self.total_genes_in_genome]) ))
            DVG.append((line_prefix+"Genes containing at least two mutants (% of all genes):", 
                        lambda summ,mutants,dataset: value_and_percentages(_N_genes_in_multiple_mutants(dataset), 
                                                                           [self.total_genes_in_genome]) ))
            DVG.append((line_prefix+"  (total genes in genome annotation data):", 
                        lambda summ,mutants,dataset: "(%s)"%self.total_genes_in_genome ))
            # MAYBE-TODO put some kind of maximum on this or collapse into ranges rather than listing all the numbers?
            for mutantN in sorted(set.union(*[set(self.get_gene_dict_by_mutant_number(dataset)) for dataset in datasets])):
                DVG.append((line_prefix+"Genes with %s mutants (%% of all genes):"%mutantN, 
                            lambda summ,mutants,dataset,N=mutantN: value_and_percentages(
                                                                     len(self.get_gene_dict_by_mutant_number(dataset)[N]), 
                                                                      [self.total_genes_in_genome]) ))
                DVG.append((line_prefix+"  (some gene names):",
                            lambda summ,mutants,dataset,N=mutantN: self._make_genelist_str(
                                                    self.get_gene_dict_by_mutant_number(dataset)[N], N_genes_to_print) ))
            # LATER-TODO Add some measure of mutations, like how many mutants have <50% perfect reads (or something - the number should probably be a command-line option).  Maybe how many mutants have <20%, 20-80%, and >80% perfect reads (or 10 and 90, or make that a variable...)
        
        ### print the data: line description and tab-separated list of values for each dataset.
        # for multi-dataset mutants, write header line with dataset names in order
        if self.multi_dataset:      
            OUTPUT.write(header_prefix + 'DATASETS\t' + '\t'.join(datasets) + '\n')
        for line_description, line_value_getter in descriptions_and_value_getters:
            OUTPUT.write(line_description)
            for summ,dataset in summaries_and_datasets:
                # TODO is this self.mutants_in_dataset thing actually needed?...
                mutants = list(self.mutants_in_dataset(dataset))
                OUTPUT.write('\t' + line_value_getter(summ,mutants,dataset))
            OUTPUT.write('\n')


    def print_data(self, OUTPUT=sys.stdout, sort_data_by=None, N_sequences=None, header_line=True, header_prefix="# "):
        """ Print full data, one line per mutant: position data, gene info, read counts, optionally sequences.
        (see the file header line for exactly what all the output fields are).

        For a normal dataset, print position/gene info (7 fields), total_reads, perfect_reads, N_sequence_variants, 
         the N_sequences most common sequences and counts (alternating, taking 2*N_sequences fields), 
         and gene annotation if it exists (arbitrary number of tab-separated fields).

        For a multi-dataset, print position/gene info (7 fields), the main sequence (taken over all the datasets),
         then reads_in_<x> and perfect_in_<x> for each dataset (total and perfect reads - total 2*N_datasets fields), 
         then gene annotation if it exists.

        Data is printed to OUTPUT, which should be an open file object (stdout by default).
         Output is tab-separated, with optional header starting with "# ".  
        """
        if self.multi_dataset and N_sequences is not None and N_sequences!=1:
            raise MutantError("Only one sequence can currently be printed in print_data for multi-datasets!")
        if N_sequences is None and not self.multi_dataset:     N_sequences = 2

        ### print the header line (different for normal and multi-datasets)
        # MAYBE-TODO should printing the gene info be optional?  Maybe... 
        #   Would save space when there isn't any meaningful gene info to print anyway.
        if header_line:
            header = ['chromosome','strand','min_position','full_position', 'gene','orientation','feature']
            if not self.multi_dataset:
                header += ['total_reads','perfect_reads', 'N_sequence_variants']
                for N in range(1,N_sequences+1):    
                    header += ['read_sequence_%s'%N, 'seq_%s_count'%N]
            else:
                header += ['main_sequence']
                for dataset_name in self.dataset_order:
                    header += ['reads_in_%s'%dataset_name, 'perfect_in_%s'%dataset_name]
            header += self.gene_annotation_header
            OUTPUT.write(header_prefix + '\t'.join(header) + "\n")

        ### sort all mutants by position or readcount (for single datasets only for now), or don't sort at all
        all_mutants = iter(self)
        if sort_data_by=='position':
            data = sorted(all_mutants, key = lambda m: m.position)
            # x.position here is an Insertion_position object and has a sensible cmp function
        elif sort_data_by=='read_count':
            if self.multi_dataset:  
                raise MutantError("Sorting by readcount in print_data not implemented for multi-datasets!")
            data = sorted(all_mutants, key = lambda m: (m.total_read_count,m.perfect_read_count,m.position), reverse=True)
        else:
            data = all_mutants

        # create "empty" annotation line with the correct number of fields, for genes that weren't in annotation file
        if self.gene_annotation_header:
            missing_gene_annotation_data = ['NO GENE DATA'] + ['' for x in self.gene_annotation_header[:-1]]

        ### for each mutant, print the mutant data line (different for normal and multi-datasets)
        for mutant in data:
            mutant_data = [mutant.position.chromosome, mutant.position.strand, mutant.position.min_position, 
                           mutant.position.full_position, mutant.gene, mutant.orientation, mutant.gene_feature] 
            if not self.multi_dataset:
                mutant_data += [mutant.total_read_count, mutant.perfect_read_count, mutant.unique_sequence_count]
                for N in range(1,N_sequences+1):
                    mutant_data += list(mutant.get_main_sequence(N))
                    # MAYBE-TODO also give the length and number of mutations for each sequence? Optionally?  
                    #   Length is easy, but do I even keep track of mutation number?  I probably should...
            else:
                mutant_data += [mutant.get_main_sequence(1)[0]]
                for dataset_name in self.dataset_order:
                    mutant_data += [mutant.by_dataset[dataset_name].total_read_count, 
                                    mutant.by_dataset[dataset_name].perfect_read_count]
            # add gene annotation, or a line with the right number of fields if gene annotation is missing
            if self.gene_annotation_header:
                if mutant.gene_annotation:  mutant_data += mutant.gene_annotation
                else:                       mutant_data += missing_gene_annotation_data
            OUTPUT.write('\t'.join([str(x) for x in mutant_data]))
            OUTPUT.write('\n')


############################################### Unit-tests ##############################################################

class Testing_position_functionality(unittest.TestCase):
    """ Unit-tests for position-related classes and functions. """

    from deepseq_utilities import Fake_deepseq_objects
    Fake_HTSeq_genomic_pos = Fake_deepseq_objects.Fake_HTSeq_genomic_pos

    def test__Insertion_position(self):
        # different ways of making the same position or approximately same position, 
        #  by defining position_after, position_before, or both
        a_pos1 = Insertion_position('chr','+', '?-101')
        a_pos2 = Insertion_position('chr','+', full_position='?-101')
        a_pos3 = Insertion_position('chr','+', position_after=101)
        a_pos4 = Insertion_position('chr','+', position_after='101')
        a_pos5 = Insertion_position('chr','+', '-101')
        all_after_positions = [a_pos1, a_pos2, a_pos3, a_pos4, a_pos5]
        b_pos1 = Insertion_position('chr','+', '100-?')
        b_pos2 = Insertion_position('chr','+', full_position='100-?')
        b_pos3 = Insertion_position('chr','+', position_before=100)
        b_pos4 = Insertion_position('chr','+', position_before='100')
        b_pos5 = Insertion_position('chr','+', '100-')
        all_before_positions = [b_pos1, b_pos2, b_pos3, b_pos4, b_pos5]
        c_pos1 = Insertion_position('chr','+', '100-101')
        c_pos2 = Insertion_position('chr','+', full_position='100-101')
        c_pos3 = Insertion_position('chr','+', position_before=100, position_after=101)
        c_pos4 = Insertion_position('chr','+', position_before='100', position_after=101)
        c_pos5 = Insertion_position('chr','+', '100-101')
        all_both_positions = [c_pos1, c_pos2, c_pos3, c_pos4, c_pos5]
        all_positions = all_after_positions + all_before_positions + all_both_positions
        # most things are exactly the same for all these positions
        assert all([pos.chromosome == 'chr' for pos in all_positions])
        assert all([pos.strand == '+' for pos in all_positions])
        assert all([pos.min_position == 100 for pos in all_positions])
        assert all([pos.max_position == 101 for pos in all_positions])
        # position_before and position_after is defined for some and not others
        assert all([pos.full_position == '?-101' for pos in all_after_positions])
        assert all([pos.full_position == '100-?' for pos in all_before_positions])
        assert all([pos.full_position == '100-101' for pos in all_both_positions])
        assert all([pos.position_before is None for pos in all_after_positions])
        assert all([pos.position_before == 100 for pos in all_before_positions+all_both_positions])
        assert all([pos.position_after is None for pos in all_before_positions])
        assert all([pos.position_after == 101 for pos in all_after_positions+all_both_positions])
        # printing - str gives just the basic info, repr gives the function to generate the object
        assert all([str(pos) == 'chr + 100-?' for pos in all_before_positions])
        assert all([repr(pos) == "Insertion_position('chr', '+', full_position='100-?', immutable=False)" 
                    for pos in all_before_positions])
        # comparison - positions are only equal if they're exactly identical, so ?-101 and 100-101 aren't equal
        assert a_pos1 == a_pos2 == a_pos3 == a_pos4 == a_pos5
        assert b_pos1 == b_pos2 == b_pos3 == b_pos4 == b_pos5
        assert c_pos1 == c_pos2 == c_pos3 == c_pos4 == c_pos5
        assert a_pos1 != b_pos1 != c_pos1 != a_pos1     # need a_pos1 twice to get it compared to both of the others
        # sorting - based on chromosome names/numbers, then min_pos, then strand
        pos1 = Insertion_position('chr','+', '10-201')
        pos2 = Insertion_position('chr','+', '?-101')
        pos3 = Insertion_position('chr','+', '200-?')
        pos4 = Insertion_position('chr','-', '200-?')
        pos5 = Insertion_position('chr2','+', '?-101')
        pos6 = Insertion_position('chr12','+', '?-101')
        pos7 = Insertion_position('other_chr','+', '?-101')
        pos8 = Insertion_position('other_chr_4','+', '?-101')
        pos9 = Insertion_position('other_chr_13','+', '?-101')
        assert pos1 < pos2 < pos3 < pos4 < pos5 < pos6 < pos7 < pos8 < pos9
        ### copying position - same contents, different ID
        pos5 = pos1.copy()
        assert pos5 == pos1
        assert pos5 is not pos1
        ### invalid creation options
        # at least one position arg is required, and there must be at least one meaningful number overall
        self.assertRaises(ValueError, Insertion_position, 'chr','+')
        self.assertRaises(ValueError, Insertion_position, 'chr','+', full_position=None, position_before=None)
        self.assertRaises(ValueError, Insertion_position, 'chr','+', position_before=None, position_after=None)
        self.assertRaises(ValueError, Insertion_position, 'chr','+', full_position='?-?')
        self.assertRaises(ValueError, Insertion_position, 'chr','+', full_position='-')
        # if providing full_position, can't provide another position arg
        self.assertRaises(ValueError, Insertion_position, 'chr','+','100-?',position_before=100)
        self.assertRaises(ValueError, Insertion_position, 'chr','+','100-?',position_after=100)
        self.assertRaises(ValueError, Insertion_position, 'chr','+','100-?',position_before=100,position_after=100)
        # full_position arg must be proper format
        for bad_fullpos in ['100', '100?', 100, (100,200), [100,200], True, False, '1-2-3', 'adaf', 'as-1']:
            self.assertRaises(ValueError, Insertion_position, 'chr','+', bad_fullpos)
        ### testing immutability/mutability and hashability
        # object is mutable and unhashable to start with unless specified otherwise
        pos1 = Insertion_position('chr','+', '0-?')
        assert str(pos1) == 'chr + 0-?'
        pos1.position_after = 5
        assert str(pos1) == 'chr + 0-5'
        pos1.new_attribute = 'test'
        assert pos1.new_attribute == 'test'
        self.assertRaises(MutantError, set, [pos1])
        self.assertRaises(MutantError, dict, [(pos1, 1)])
        # if we make it immutable, or initialize it as immutable from the start, it's immutable and hashable
        pos1.make_immutable()
        pos2 = Insertion_position('chr','+', '0-?', immutable=True)
        for pos in (pos1,pos2):
            # have to use execute workaround here to use self.assertRaises on statements instead of expressions;
            #  if this was python 2.7 I could just do "with self.assertRaises(MutantError): <bunch of statements>"
            def execute(S, context):  exec S in context
            self.assertRaises(MutantError, execute, "pos.chromosome = 'chr3'", locals())
            self.assertRaises(MutantError, execute, "pos.strand = '-'", locals())
            self.assertRaises(MutantError, execute, "pos.position_after = '5'", locals())
            self.assertRaises(MutantError, execute, "pos.position_before = '4'", locals())
            self.assertRaises(MutantError, execute, "pos.new_attribute = '4'", locals())
            set([pos1])
            dict([(pos1, 1)])
        # if we make it mutable, it's mutable and unhashable again
        pos1.make_mutable_REMEMBER_CLEANUP_FIRST()
        pos2.make_mutable_REMEMBER_CLEANUP_FIRST()
        for pos in pos1,pos2:
            pos.position_before = None
            pos.position_after = 500
            assert str(pos) == 'chr + ?-500'
            pos.new_attribute = 'test'
            assert pos.new_attribute == 'test'
            self.assertRaises(MutantError, set, [pos])
            self.assertRaises(MutantError, dict, [(pos, 1)])


    def test__get_insertion_pos_from_HTSeq_position(self):
        # should raise exception for invalid HTSeq_pos argument
        for bad_HTSeq_pos in [None, '', 'aaa', 0, 1, 0.65, [], {}, True, False]:
            for cassette_end in SEQ_ENDS:
                self.assertRaises(ValueError, get_insertion_pos_from_HTSeq_read_pos, None, cassette_end)
        # should raise exception for invalid cassette_end
        fake_pos = self.Fake_HTSeq_genomic_pos('C', '+', 0, 5)
        for bad_cassette_end in ['','aaa',0,1,[],{},None,True,False,'start','end','middle','read','leftmost','rightmost']:
            self.assertRaises(ValueError, get_insertion_pos_from_HTSeq_read_pos, fake_pos, bad_cassette_end)
        ### testing normal functionality: should return an Insertion_position instance with the same chromosome, 
        #    and strand/position depending on the arguments in a somewhat complicated way.
        ## a few spot-checks outside the loop 
        #   (based on the example at the end of "Possible read sides/directions" section in ../notes.txt)
        # remember HTSeq position is 0-based and end-exclusive, and the position I want is 1-based end-inclusive!  
        #  So in the end the two relevant 1-based numbers end up being the same as the 0-based positions,
        #  because in the case of the start, min_position is actually start-1, and in the case of the end, we're adding 1
        #  to switch from 0-based to 1-based but then subtracting 1 to make the end the last base instead of the base after
        fake_pos_plus = self.Fake_HTSeq_genomic_pos('C', '+', 3, 7)
        fake_pos_minus = self.Fake_HTSeq_genomic_pos('C', '-', 3, 7)
        assert get_insertion_pos_from_HTSeq_read_pos(fake_pos_plus, '5prime',reads_are_reverse=False).min_position == 7
        assert get_insertion_pos_from_HTSeq_read_pos(fake_pos_plus, '3prime',reads_are_reverse=False).min_position == 3
        assert get_insertion_pos_from_HTSeq_read_pos(fake_pos_minus,'5prime',reads_are_reverse=True ).min_position == 7
        assert get_insertion_pos_from_HTSeq_read_pos(fake_pos_minus,'3prime',reads_are_reverse=True ).min_position == 3
        assert get_insertion_pos_from_HTSeq_read_pos(fake_pos_plus, '5prime',reads_are_reverse=True ).min_position == 3
        assert get_insertion_pos_from_HTSeq_read_pos(fake_pos_plus, '3prime',reads_are_reverse=True ).min_position == 7
        assert get_insertion_pos_from_HTSeq_read_pos(fake_pos_minus,'5prime',reads_are_reverse=False).min_position == 3
        assert get_insertion_pos_from_HTSeq_read_pos(fake_pos_minus,'3prime',reads_are_reverse=False).min_position == 7
        # more checks in a loop - see above for what we expect!
        for (strand_in, if_reverse, strand_out) in [('+',False,'+'), ('-',False,'-'), ('+',True,'-'), ('-',True,'+')]:
            for (start,end) in [(0,5), (0,100), (10,11), (5,44)]:
                fake_pos = self.Fake_HTSeq_genomic_pos('C', strand_in, start, end)
                result_5prime = get_insertion_pos_from_HTSeq_read_pos(fake_pos, '5prime', if_reverse)
                result_3prime = get_insertion_pos_from_HTSeq_read_pos(fake_pos, '3prime', if_reverse)
                assert result_5prime.chromosome == result_3prime.chromosome == 'C'
                assert result_5prime.strand == result_3prime.strand == strand_out
                if strand_out=='+':
                    assert result_5prime.min_position == end
                    assert result_3prime.min_position == start
                elif strand_out=='-':
                    assert result_5prime.min_position == start
                    assert result_3prime.min_position == end


class Testing_Insertional_mutant(unittest.TestCase):
    """ Unit-tests for the Insertional_mutant class and its methods. """

    def test__init(self):
        for chromosome in ['chr1', 'chromosome_2', 'chrom3', 'a', 'adfads', '100', 'scaffold_88']:
            for strand in ['+','-']:
                for position in [1,2,5,100,10000,4323423]:
                    ins_pos_5prime = Insertion_position(chromosome,strand,position_before=position)
                    ins_pos_3prime = Insertion_position(chromosome,strand,position_after=position)
                    # test "normal" mutants - check all details, including position
                    mutant_5prime = Insertional_mutant(ins_pos_5prime)
                    mutant_3prime = Insertional_mutant(ins_pos_3prime)
                    mutant_readcount_only = Insertional_mutant(ins_pos_5prime, readcount_related_only=True)
                    mutant_multi_dataset = Insertional_mutant(ins_pos_5prime, multi_dataset=True)
                    # test position details (only for the two "normal" mutants)
                    assert mutant_5prime.position.min_position == position
                    assert mutant_3prime.position.min_position == position-1
                    assert mutant_5prime.position.max_position == position+1
                    assert mutant_3prime.position.max_position == position
                    assert mutant_5prime.position.full_position == "%s-?"%(position)
                    assert mutant_3prime.position.full_position == "?-%s"%position
                    # test non-readcount-related info for all mutants except mutant_readcount_only
                    for mutant in [mutant_5prime, mutant_3prime, mutant_multi_dataset]:
                        assert mutant.position.chromosome == chromosome
                        assert mutant.position.strand == strand
                        assert mutant.gene == SPECIAL_GENE_CODES.not_determined
                        assert mutant.orientation == '?'
                        assert mutant.gene_feature == '?'
                    # test readcount-related info for all mutants except mutant_multi_dataset
                    for mutant in [mutant_5prime, mutant_3prime, mutant_readcount_only]:
                        assert mutant.total_read_count == 0
                        assert mutant.perfect_read_count == 0
                        assert mutant.unique_sequence_count == 0
                        assert mutant.sequences_and_counts == {}
                    # test readcount-related info for mutant_multi_dataset
                    assert all([x.total_read_count == 0 for x in mutant_multi_dataset.by_dataset.values()])
                    assert all([x.perfect_read_count == 0 for x in mutant_multi_dataset.by_dataset.values()])
                    assert all([x.unique_sequence_count == 0 for x in mutant_multi_dataset.by_dataset.values()])
                    assert all([x.sequences_and_counts == {} for x in mutant_multi_dataset.by_dataset.values()])

    def test__add_read(self):
        # using fake HTSeq alignment class from deepseq_utilities; defining one perfect and one imperfect alignment
        # note: the detailed mutation-counting methods are imported from deepseq_utilities and unit-tested there.
        from deepseq_utilities import Fake_deepseq_objects
        Fake_HTSeq_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment
        perfect_aln = Fake_HTSeq_alignment(seq='AAA', cigar_string='===')  # CIGAR = means match
        imperfect_aln = Fake_HTSeq_alignment(seq='GGG', cigar_string='=X=')  # CIGAR X means mismatch
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        # adding perfect and imperfect to mutant increases all the counts as expected
        mutant.add_read(perfect_aln, read_count=3)
        assert mutant.total_read_count == mutant.perfect_read_count == 3
        assert mutant.unique_sequence_count == 1
        assert mutant.sequences_and_counts == {'AAA':3}
        mutant.add_read(imperfect_aln, read_count=1)
        assert mutant.total_read_count == 4
        assert mutant.perfect_read_count == 3
        assert mutant.unique_sequence_count == 2
        assert mutant.sequences_and_counts == {'AAA':3, 'GGG':1}
        # it should be impossible to add a read to a specific dataset in a single-dataset mutant 
        self.assertRaises(MutantError, mutant.add_read, perfect_aln, read_count=3, dataset_name='d1')
        # same for a multi-dataset mutant - this time we need to specify which dataset we're adding to
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3), multi_dataset=True)
        assert len(mutant.by_dataset) == 0
        mutant.add_read(perfect_aln, read_count=3, dataset_name='d1')
        assert len(mutant.by_dataset) == 1
        assert mutant.by_dataset['d1'].total_read_count == mutant.by_dataset['d1'].perfect_read_count == 3
        assert mutant.by_dataset['d1'].unique_sequence_count == 1
        assert mutant.by_dataset['d1'].sequences_and_counts == {'AAA':3}
        mutant.add_read(imperfect_aln, read_count=1, dataset_name='d1')
        assert len(mutant.by_dataset) == 1
        assert mutant.by_dataset['d1'].total_read_count == 4
        assert mutant.by_dataset['d1'].perfect_read_count == 3
        assert mutant.by_dataset['d1'].unique_sequence_count == 2
        assert mutant.by_dataset['d1'].sequences_and_counts == {'AAA':3, 'GGG':1}
        # now adding a read to another dataset - nothing changes in dataset d1, but we have new dataset d2 numbers
        mutant.add_read(imperfect_aln, read_count=1, dataset_name='d2')
        assert len(mutant.by_dataset) == 2
        assert mutant.by_dataset['d1'].total_read_count == 4
        assert mutant.by_dataset['d2'].total_read_count == 1
        assert mutant.by_dataset['d2'].perfect_read_count == 0
        assert mutant.by_dataset['d2'].unique_sequence_count == 1
        assert mutant.by_dataset['d2'].sequences_and_counts == {'GGG':1}
        # it should be impossible to add a read to a multi-dataset mutant without giving a dataset_name
        self.assertRaises(MutantError, mutant.add_read, perfect_aln, read_count=3)

    def test__update_gene_info(self):
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        assert mutant.gene == SPECIAL_GENE_CODES.not_determined
        assert mutant.orientation == mutant.gene_feature == '?'
        # updating no-info mutant with no info - no change
        mutant.update_gene_info(SPECIAL_GENE_CODES.not_determined, '?', '?')
        assert mutant.gene == SPECIAL_GENE_CODES.not_determined
        assert mutant.orientation == mutant.gene_feature == '?'
        # updating no-info mutant with useful info - update goes through
        mutant.update_gene_info('gene1', '+', 'f')
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        # updating info-containing mutant with no info - no change
        mutant.update_gene_info(SPECIAL_GENE_CODES.not_determined, '?', '?')
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        # updating info-containing mutant with same info - no change
        mutant.update_gene_info('gene1', '+', 'f')
        assert mutant.gene == 'gene1'
        assert mutant.orientation == '+'
        assert mutant.gene_feature == 'f'
        # updating info-containig mutant with OTHER info - exception
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene2', '+', 'f')
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene1', '-', 'f')
        self.assertRaises(MutantError, mutant.update_gene_info, 'gene1', '+', 'g')


    def test__merge_mutant(self):
        mutant1 = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant1.add_counts(2,2,1)
        mutant1.add_sequence_and_counts('AAA',2)
        mutant2 = Insertional_mutant(Insertion_position('chr','+',position_before=5))
        mutant2.add_counts(1,0,1)
        mutant2.add_sequence_and_counts('AA',1)
        # merging two mutants - read counts act as expected
        assert (mutant1.total_read_count, mutant1.perfect_read_count) == (2, 2)
        assert (mutant2.total_read_count, mutant2.perfect_read_count) == (1, 0)
        assert (mutant1.sequences_and_counts, mutant2.sequences_and_counts) == ({'AAA':2}, {'AA':1})
        mutant1.merge_mutant(mutant2)
        assert (mutant1.total_read_count, mutant1.perfect_read_count) == (3, 2)
        assert (mutant2.total_read_count, mutant2.perfect_read_count) == (0, 0)
        assert (mutant1.sequences_and_counts, mutant2.sequences_and_counts) == ({'AAA':2, 'AA':1}, {})
        # if we define the gene for one of the mutants but not the other, still works
        #  (note that starting here mutant2 has 0 reads, but we can still merge it all we want just to check for errors)
        mutant2.gene = 'X'
        mutant1.merge_mutant(mutant2, check_gene_data=True)
        # if we define the gene for both mutants and it's the same gene, also works
        mutant1.gene = 'X'
        mutant1.merge_mutant(mutant2, check_gene_data=True)
        # if we define the gene for both mutants and it's a DIFFERENT gene, we get an exception unless we don't check
        mutant2.gene = 'Y'
        self.assertRaises(MutantError, mutant1.merge_mutant, mutant2, check_gene_data=True)
        mutant1.merge_mutant(mutant2, check_gene_data=False)
        # mutant-merging for multi-dataset mutants currently NOT IMPLEMENTED

    def test__add_counts(self):
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant.add_counts(0,0,0)
        assert mutant.total_read_count == 0
        assert mutant.perfect_read_count == 0
        assert mutant.sequences_and_counts == {}
        mutant.add_counts(2,2,1)
        assert mutant.total_read_count == 2
        assert mutant.perfect_read_count == 2
        assert mutant.sequences_and_counts == {}
        mutant.add_counts(2,1,1)
        assert mutant.total_read_count == 4
        assert mutant.perfect_read_count == 3
        assert mutant.sequences_and_counts == {}
        # how mutant.unique_sequence_count changes depends on assume_new_sequences:
        #  - if False, mutant.unique_sequence_count is max(new_seq_count, mutant.unique_sequence_count)
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant.add_counts(0,0,0,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 0
        mutant.add_counts(1,1,1,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 1
        mutant.add_counts(2,2,2,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 2
        mutant.add_counts(1,1,1,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 2
        mutant.add_counts(2,2,2,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 2
        #  - if True and mutant.unique_sequence_count is new_seq_count + mutant.unique_sequence_count
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant.add_counts(0,0,0,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 0
        mutant.add_counts(1,1,1,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 1
        mutant.add_counts(2,2,2,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 3
        mutant.add_counts(1,1,1,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 4
        mutant.add_counts(2,2,2,assume_new_sequences=True)
        assert mutant.unique_sequence_count == 6
        mutant.add_counts(2,2,2,assume_new_sequences=False)
        assert mutant.unique_sequence_count == 6
        ### works for multi-dataset mutants (not looking at too much detail - should work the same as for single)
        multi_mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3), multi_dataset=True)
        multi_mutant.add_counts(2,2,1, dataset_name='d')
        assert multi_mutant.by_dataset['d'].total_read_count == 2
        assert multi_mutant.by_dataset['d'].perfect_read_count == 2
        assert multi_mutant.by_dataset['d'].unique_sequence_count == 1
        # doesn't work for multi-dataset mutants if dataset_name not given, or for single if given
        self.assertRaises(MutantError, multi_mutant.add_counts, 1,1,1)
        self.assertRaises(MutantError, mutant.add_counts, 1,1,1, dataset_name='d1')

    def test__add_sequence_and_counts(self):
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        # adding sequence/count to mutant.sequences_and_counts, WITHOUT touching mutant.unique_sequence_count
        mutant.add_sequence_and_counts('AAA',2,add_to_uniqseqcount=False)
        assert mutant.sequences_and_counts == {'AAA':2}
        assert mutant.unique_sequence_count == 0
        # adding sequence/count to mutant.sequences_and_counts, and incrementing mutant.unique_sequence_count if warranted:
        #  - if adding a sequence that was already there, don't increment
        mutant.add_sequence_and_counts('AAA',2,add_to_uniqseqcount=True)
        assert mutant.sequences_and_counts == {'AAA':4}
        assert mutant.unique_sequence_count == 0
        #  - if adding a new sequence, increment
        mutant.add_sequence_and_counts('GGG',2,add_to_uniqseqcount=True)
        assert mutant.sequences_and_counts == {'AAA':4, 'GGG':2}
        assert mutant.unique_sequence_count == 1
        #  - if adding a new sequence but mutant.unique_sequence_count is already higher than expected, don't increment
        mutant.unique_sequence_count = 5
        mutant.add_sequence_and_counts('CCC',2,add_to_uniqseqcount=True)
        assert mutant.sequences_and_counts == {'AAA':4, 'GGG':2, 'CCC':2}
        assert mutant.unique_sequence_count == 5
        # make sure it raises an error if given a non-numeric argument
        for not_a_number in [None,'','a','GGG']:
            self.assertRaises(TypeError,mutant.add_sequence_and_counts,'CCC',not_a_number)
        ### works for multi-dataset mutants (not looking at too much detail - should work the same as for single)
        multi_mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3), multi_dataset=True)
        assert multi_mutant.by_dataset['d'].sequences_and_counts == {}
        multi_mutant.add_sequence_and_counts('GGG',2, dataset_name='d')
        assert multi_mutant.by_dataset['d'].sequences_and_counts == {'GGG':2}
        # doesn't work for multi-dataset mutants if dataset_name not given, or for single if given
        self.assertRaises(MutantError, multi_mutant.add_sequence_and_counts, 'GGG',1)
        self.assertRaises(MutantError, mutant.add_sequence_and_counts, 'GGG',1, dataset_name='d1')

    def test__get_main_sequence(self):
        # single-dataset mutant
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        assert mutant.get_main_sequence() == ('',0)
        assert mutant.get_main_sequence(1) == ('',0)
        assert mutant.get_main_sequence(4) == ('',0)
        mutant.add_sequence_and_counts('AAA',1)
        mutant.add_sequence_and_counts('GGG',2)
        assert mutant.sequences_and_counts == {'AAA':1, 'GGG':2}
        assert mutant.get_main_sequence() == ('GGG',2)
        assert mutant.get_main_sequence(1) == ('GGG',2)
        assert mutant.get_main_sequence(2) == ('AAA',1)
        assert mutant.get_main_sequence(3) == ('',0)
        assert mutant.get_main_sequence(4) == ('',0)
        mutant.add_sequence_and_counts('CCC',1)
        mutant.add_sequence_and_counts('AAA',2)
        assert mutant.sequences_and_counts == {'AAA':3, 'GGG':2, 'CCC':1}
        assert mutant.get_main_sequence() == ('AAA',3)
        assert mutant.get_main_sequence(1) == ('AAA',3)
        assert mutant.get_main_sequence(2) == ('GGG',2)
        assert mutant.get_main_sequence(3) == ('CCC',1)
        assert mutant.get_main_sequence(4) == ('',0)
        assert mutant.get_main_sequence(5) == ('',0)
        # multi-dataset mutant - getting the top sequence from single dataset or from all datasets added together
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3), multi_dataset=True)
        mutant.add_sequence_and_counts('AAA',3, dataset_name='d1')
        mutant.add_sequence_and_counts('GGG',2, dataset_name='d1')
        mutant.add_sequence_and_counts('TTT',3, dataset_name='d2')
        mutant.add_sequence_and_counts('GGG',2, dataset_name='d2')
        assert mutant.get_main_sequence(1, dataset_name='d1') == ('AAA',3)
        assert mutant.get_main_sequence(1, dataset_name='d2') == ('TTT',3)
        assert mutant.get_main_sequence(1) == ('GGG',4)     # GGG is the most common sequence if we add both datasets

    def test__convert_to_multi_dataset(self):
        # converting to multi-dataset WITH saving current data to dataset d
        mutant1 = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant1.add_counts(2,1,1)
        mutant1.add_sequence_and_counts('AAA',2)
        mutant1.convert_to_multi_dataset(current_dataset_name='d')
        assert len(mutant1.by_dataset) == 1
        assert mutant1.by_dataset['d'].total_read_count == 2
        assert mutant1.by_dataset['d'].perfect_read_count == 1
        assert mutant1.by_dataset['d'].unique_sequence_count == 1
        assert mutant1.by_dataset['d'].sequences_and_counts == {'AAA':2}
        # converting to multi-dataset WITHOUT saving current data - the values for any dataset come out as 0/empty 
        mutant2 = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant2.add_counts(2,1,1)
        mutant2.add_sequence_and_counts('AAA',2)
        mutant2.convert_to_multi_dataset()
        assert all([dataset.total_read_count == 0 for dataset in mutant2.by_dataset.values()])
        assert all([dataset.perfect_read_count == 0 for dataset in mutant2.by_dataset.values()])
        assert all([dataset.unique_sequence_count == 0 for dataset in mutant2.by_dataset.values()])
        assert all([dataset.sequences_and_counts == {} for dataset in mutant2.by_dataset.values()])
        # converting an already multi-dataset mutant to multi-dataset: MutantError if not ignore_if_already_multi, 
        #  otherwise works and doesn't change anything - all the same asserts
        self.assertRaises(MutantError, mutant1.convert_to_multi_dataset, current_dataset_name='d')
        self.assertRaises(MutantError, mutant2.convert_to_multi_dataset)
        mutant1.convert_to_multi_dataset(current_dataset_name='d', ignore_if_already_multi=True)
        assert len(mutant1.by_dataset) == 1
        assert mutant1.by_dataset['d'].total_read_count == 2
        mutant2.convert_to_multi_dataset(ignore_if_already_multi=True)
        assert mutant2.by_dataset['d'].total_read_count == 0
        assert all([dataset.total_read_count == 0 for dataset in mutant2.by_dataset.values()])

    def test__add_other_mutant_as_dataset(self):
        mutant1 = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant1.add_counts(2,1,1)
        mutant1.add_sequence_and_counts('AAA',2)
        mutant2 = Insertional_mutant(Insertion_position('chr','+',position_before=3))
        mutant2.add_counts(1,0,1)
        mutant2.add_sequence_and_counts('GGG',1)
        # adding a mutant to a single-dataset mutant should fail
        self.assertRaises(MutantError, mutant1.add_other_mutant_as_dataset, mutant2, 'd2')
        # adding a mutant to a multi-dataset mutant should work
        mutant1.convert_to_multi_dataset(current_dataset_name='d1')
        mutant1.add_other_mutant_as_dataset(mutant2, 'd2')
        assert mutant1.by_dataset['d1'].total_read_count == 2
        assert mutant1.by_dataset['d1'].perfect_read_count == 1
        assert mutant1.by_dataset['d1'].unique_sequence_count == 1
        assert mutant1.by_dataset['d1'].sequences_and_counts == {'AAA':2}
        assert mutant1.by_dataset['d2'].total_read_count == 1
        assert mutant1.by_dataset['d2'].perfect_read_count == 0
        assert mutant1.by_dataset['d2'].unique_sequence_count == 1
        assert mutant1.by_dataset['d2'].sequences_and_counts == {'GGG':1}
        # can't add overwrite an existing dataset name, unless overwrite==True
        self.assertRaises(MutantError, mutant1.add_other_mutant_as_dataset, mutant2, 'd2')
        mutant1.add_other_mutant_as_dataset(mutant2, 'd2', overwrite=True)
        # if the two mutants have different positions, it should fail, unless check_constant_data=False
        mutant3 = Insertional_mutant(Insertion_position('chr','+',position_before=5))
        mutant4 = Insertional_mutant(Insertion_position('chr2','+',position_before=3))
        mutant5 = Insertional_mutant(Insertion_position('chr','-',position_before=3))
        self.assertRaises(MutantError, mutant1.add_other_mutant_as_dataset, mutant3, 'd3', check_constant_data=True)
        self.assertRaises(MutantError, mutant1.add_other_mutant_as_dataset, mutant4, 'd4', check_constant_data=True)
        self.assertRaises(MutantError, mutant1.add_other_mutant_as_dataset, mutant5, 'd5', check_constant_data=True)
        mutant1.add_other_mutant_as_dataset(mutant3, 'd3', check_constant_data=False)
        mutant1.add_other_mutant_as_dataset(mutant4, 'd4', check_constant_data=False)
        mutant1.add_other_mutant_as_dataset(mutant5, 'd5', check_constant_data=False)

    def test__give_single_dataset_mutant(self):
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3), multi_dataset=True)
        mutant.add_counts(2,1,1, dataset_name='d1')
        mutant.add_sequence_and_counts('AAA',2, dataset_name='d1')
        mutant.add_counts(1,0,1, dataset_name='d2')
        mutant.add_sequence_and_counts('GGG',1, dataset_name='d2')
        # extracting two mutants and checking the values
        new_mutant_1 = mutant.give_single_dataset_mutant('d1')
        new_mutant_2 = mutant.give_single_dataset_mutant('d2')
        for new_mutant in (new_mutant_1, new_mutant_2):
            assert new_mutant.position == mutant.position
            assert new_mutant.gene == mutant.gene
        assert new_mutant_1.total_read_count == 2
        assert new_mutant_1.perfect_read_count == 1
        assert new_mutant_1.unique_sequence_count == 1
        assert new_mutant_1.sequences_and_counts == {'AAA':2}
        assert new_mutant_2.total_read_count == 1
        assert new_mutant_2.perfect_read_count == 0
        assert new_mutant_2.unique_sequence_count == 1
        assert new_mutant_2.sequences_and_counts == {'GGG':1}
        # trying to extract an inexistent dataset should fail, unless force==True
        self.assertRaises(MutantError, mutant.give_single_dataset_mutant, 'd0', force=False)
        new_mutant_0 = mutant.give_single_dataset_mutant('d0', force=True)
        assert new_mutant_0.position == mutant.position
        assert new_mutant_0.gene == mutant.gene
        assert new_mutant_0.total_read_count == 0
        assert new_mutant_0.perfect_read_count == 0
        assert new_mutant_0.unique_sequence_count == 0
        assert new_mutant_0.sequences_and_counts == {}
        # trying to extract a dataset from a single-dataset mutant should fail
        mutant1 = Insertional_mutant(Insertion_position('chr','+',position_before=3), multi_dataset=False)
        self.assertRaises(MutantError, mutant1.give_single_dataset_mutant, 'd2')

    def test__give_all_single_dataset_mutants(self):
        mutant = Insertional_mutant(Insertion_position('chr','+',position_before=3), multi_dataset=True)
        mutant.add_counts(2,1,1, dataset_name='d1')
        mutant.add_sequence_and_counts('AAA',2, dataset_name='d1')
        mutant.add_counts(1,0,1, dataset_name='d2')
        mutant.add_sequence_and_counts('GGG',1, dataset_name='d2')
        # extracting two mutants and checking the values
        all_single_mutants = mutant.give_all_single_dataset_mutants()
        assert len(all_single_mutants) == 2
        new_mutant_1 = all_single_mutants['d1']
        new_mutant_2 = all_single_mutants['d2']
        for new_mutant in (new_mutant_1, new_mutant_2):
            assert new_mutant.position == mutant.position
            assert new_mutant.gene == mutant.gene
        assert new_mutant_1.total_read_count == 2
        assert new_mutant_1.perfect_read_count == 1
        assert new_mutant_1.unique_sequence_count == 1
        assert new_mutant_1.sequences_and_counts == {'AAA':2}
        assert new_mutant_2.total_read_count == 1
        assert new_mutant_2.perfect_read_count == 0
        assert new_mutant_2.unique_sequence_count == 1
        assert new_mutant_2.sequences_and_counts == {'GGG':1}


class Testing_Insertional_mutant_pool_dataset(unittest.TestCase):
    """ Unit-tests for the Insertional_mutant_pool_dataset class and its methods. """

    def test__init(self):
        for cassette_end in SEQ_ENDS+['?']:
            for reads_are_reverse in [True,False,'?']:
                data = Insertional_mutant_pool_dataset(cassette_end=cassette_end, reads_are_reverse=reads_are_reverse)
                assert data.summary.cassette_end == cassette_end
                assert data.summary.reads_are_reverse == reads_are_reverse
                assert len(data) == 0
                assert data.summary.processed_read_count == data.summary.aligned_read_count\
                        == data.summary.perfect_read_count == 0
                assert data.summary.non_aligned_read_count == 0
                assert data.summary.discarded_read_count == 'unknown'
                assert data.summary.ignored_region_read_counts == {}
                assert sum(data.summary.strand_read_counts.values()) == 0
                assert data.summary.mutants_in_genes == data.summary.mutants_not_in_genes\
                        == data.summary.mutants_undetermined == 0
                assert data.summary.mutant_counts_by_orientation == data.summary.mutant_counts_by_feature == {}
        for cassette_end in [True, False, None, 0, 1, 0.11, 23, 'asdfas', '', 'something', [2,1], {}]:
            self.assertRaises(ValueError, Insertional_mutant_pool_dataset, cassette_end, '?')
        for reads_are_reverse in ['forward', 'reverse', None, 0.11, 23, 'asdfas', '', 'something', [2,1], {}]:
            # note that this list doesn't include 0/1 because 0==False and 1==True, and that's what "in" tests 
            self.assertRaises(ValueError, Insertional_mutant_pool_dataset, '?', reads_are_reverse)

    # LATER-TODO add unit-test for add_discarded_reads, find_genes_for_mutants, most_common_mutants, 

    def test__add_alignment_reader_to_data(self):
        pass
        # MAYBE-TODO implement using a mock-up of HTSeq_alignment?  (see Testing_single_functions for how I did that)
        #   make sure it fails if self.cassette_end isn't defined...

    def test__print_summary(self):
        pass
        # MAYBE-TODO implement based on stuff in test_data, like do_test_run in deepseq_count_alignments.py?

    def test__print_data(self):
        pass
        # MAYBE-TODO implement based on stuff in test_data, like do_test_run in deepseq_count_alignments.py?

    def _check_infile1(self, data):
        """ Check that data is as expected in test_data/count-aln__cassette-end-5prime.txt. """
        # summary
        assert data.summary.processed_read_count == 40
        assert data.summary.aligned_read_count == 30
        assert data.summary.non_aligned_read_count == 10
        assert data.summary.perfect_read_count == 22
        assert data.summary.strand_read_counts == {'+':27, '-':3}
        assert len(data) == 17
        assert data.summary.mutants_in_genes == data.summary.mutants_not_in_genes == 0
        assert data.summary.mutant_counts_by_orientation == {}
        assert data.summary.mutant_counts_by_feature == {}
        for mutant in data:
            assert mutant.gene == SPECIAL_GENE_CODES.not_determined
            assert mutant.orientation == mutant.gene_feature == '?'
        # just spot-checking some of the outputs
        mutant = data.get_mutant('reads_2_seqs_1','+',position_before=204)
        assert mutant.position.chromosome == 'reads_2_seqs_1'
        assert mutant.position.min_position == 204
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 2
        assert mutant.perfect_read_count == 2 
        assert mutant.unique_sequence_count == 1
        mutant = data.get_mutant('mutation_yes','+',position_before=204)
        assert mutant.position.chromosome == 'mutation_yes'
        assert mutant.position.min_position == 204
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 6
        assert mutant.perfect_read_count == 0
        assert mutant.unique_sequence_count == 1
        mutant = data.get_mutant('strandedness_-_normal_+_reverse','-',position_after=101)
        assert mutant.position.chromosome == 'strandedness_-_normal_+_reverse'
        assert mutant.position.min_position == 100
        assert mutant.position.strand == '-'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1
        assert mutant.unique_sequence_count == 1

    def _check_infile2(self, data2):
        """ Check that data is as expected in test_data/count-aln__cassette-end-5prime.txt. """
        # summary
        assert data2.summary.processed_read_count == data2.summary.aligned_read_count\
                == data2.summary.perfect_read_count == 40
        assert data2.summary.non_aligned_read_count == 0
        assert data2.summary.strand_read_counts == {'+':38, '-':2}
        assert len(data2) == 40
        assert data2.summary.mutants_in_genes == 39
        assert data2.summary.mutants_not_in_genes == 1
        assert data2.summary.mutants_undetermined == 0
        assert data2.summary.mutant_counts_by_orientation == {'sense':37, 'antisense':2}
        assert data2.summary.mutant_counts_by_feature['CDS'] == 6
        assert data2.summary.mutant_counts_by_feature['??'] == 2
        assert data2.summary.mutant_counts_by_feature['CDS/three_prime_UTR'] == 1
        # mutant spot-checks
        mutant = data2.get_mutant('chromosome_A','+',position_before=20)
        assert mutant.position.chromosome == 'chromosome_A'
        assert mutant.position.min_position == 20
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1 
        assert mutant.unique_sequence_count == 1
        assert mutant.gene == SPECIAL_GENE_CODES.not_found
        assert mutant.orientation == mutant.gene_feature == '-'
        mutant = data2.get_mutant('chromosome_A','+',position_before=150)
        assert mutant.position.chromosome == 'chromosome_A'
        assert mutant.position.min_position == 150
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1 
        assert mutant.unique_sequence_count == 1
        assert mutant.gene == "test.geneA0_proper_plus"
        assert mutant.orientation == 'sense'
        assert mutant.gene_feature == 'five_prime_UTR'

    def test__read_data_from_file(self):
        """ uses _check_infile1 and _check_infile2 functions for detail. """
        ## 1. input file with no gene information but more variation in other features
        input_file = 'test_data/count-aln__cassette-end-5prime.txt'
        data = Insertional_mutant_pool_dataset()
        data.read_data_from_file(input_file)
        self._check_infile1(data)
        ## 2. input file with gene information
        input_file2 = 'test_data/count-aln__with-gene-info_merged.txt'
        data2 = Insertional_mutant_pool_dataset()
        data2.read_data_from_file(input_file2)
        self._check_infile2(data2)
        ## 3. adding more data to a file that already has some...
        data.read_data_from_file(input_file, assume_new_sequences=False)
        mutant = data.get_mutant('reads_2_seqs_1','+',position_before=204)
        assert mutant.position.chromosome == 'reads_2_seqs_1'
        assert mutant.position.min_position == 204
        assert mutant.total_read_count == 4
        assert mutant.perfect_read_count == 4
        # how mutant.unique_sequence_count should act in this case depends on the value of assume_new_sequences
        assert mutant.unique_sequence_count == 1
        data.read_data_from_file(input_file, assume_new_sequences=True)
        assert mutant.unique_sequence_count == 2

    def test__pickle_unpickle(self):
        """ uses _check_infile1 and _check_infile2 functions for detail. """
        import pickle, os
        picklefile = 'test.pickle'
        ## 1. input file with no gene information but more variation in other features
        input_file = 'test_data/count-aln__cassette-end-5prime.txt'
        data = Insertional_mutant_pool_dataset(infile=input_file)
        self._check_infile1(data)
        with open(picklefile, 'w') as PICKLE:    pickle.dump(data, PICKLE)
        with open(picklefile) as PICKLE:         data_new = pickle.load(PICKLE)
        self._check_infile1(data_new)
        os.unlink(picklefile)
        ## 2. input file with gene information
        input_file2 = 'test_data/count-aln__with-gene-info_merged.txt'
        data2 = Insertional_mutant_pool_dataset(infile=input_file2)
        self._check_infile2(data2)
        with open(picklefile, 'w') as PICKLE:    pickle.dump(data2, PICKLE)
        with open(picklefile) as PICKLE:         data2_new = pickle.load(PICKLE)
        self._check_infile2(data2_new)
        os.unlink(picklefile)


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    # if I wanted more control I could do this instead:
    #import os
    #unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromName(os.path.splitext(sys.argv[0])[0]))
    #   (autodetection of all tests - see http://docs.python.org/library/unittest.html#unittest.TestLoader)
    # there's probably also some way to easily get all tests from the current file without passing the name, but I haven't found it yet...
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main(argv=[sys.argv[0]])

