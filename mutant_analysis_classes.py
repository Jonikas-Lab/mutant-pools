#!/usr/bin/env python
"""
Module containing classes and functions for analysis of deepseq data related to insertional mutant libraries.

This is a module to be imported and used by other programs.  Running it directly runs the built-in unit-test suite.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011
"""

# basic libraries
import sys, re
import unittest
from collections import defaultdict
from itertools import combinations
# other libraries
import HTSeq
from BCBio import GFF
from Bio import SeqFeature
# my modules
from general_utilities import split_into_N_sets_by_counts, add_dicts_of_ints
from DNA_basic_utilities import SEQ_ENDS, SEQ_STRANDS, SEQ_DIRECTIONS, SEQ_ORIENTATIONS, position_test_contains, position_test_overlap
from deepseq_utilities import get_seq_count_from_collapsed_header, check_mutation_count_try_all_methods


class SPECIAL_GENE_CODES(object):
    not_determined = "gene_unknown"
    chromosome_not_in_reference = "unknown_chrom"
    not_found = "no_gene_found"
# it seems like I have to set SPECIAL_GENE_CODES.all_codes afterward because I can't access __dict__ from inside the class
SPECIAL_GENE_CODES.all_codes = [value for (name,value) in SPECIAL_GENE_CODES.__dict__.items() if not name.startswith('__')]


### Functions/classes for dealing with alignment/insertion positions

# has to be a new-style object-based class due to the immutability/hashability thing
class Insertion_position(object):
    """ A descriptor of the position of a genomic insertion, with handling of before/after sides and ambiguity.

    Attributes: 
     - chromosome and strand - the chromosome the insertion is on, and the strand it's in sense orientation to.
     - position_before and position_after - positions before and after the insertion site (1-based): 
         integers, or None if unknown;
     - min_position and max_position - lowest/highest possible position values as plain numbers, no ambiguity
     - full_position - a string describing the full position: 3-4 for exact positions, 3-? or 4-? if one side is unknown. 

    Methods: 
     - comparison/sorting: __cmp__ is based on chromosome name/number, min_/max_position, strand, position_before/_after, 
                            in that order - see __cmp__ docstring for more details/explanations
     - copy method returns an identical but separate copy of the object (not just another reference to the same object)
     - printing: __str__ returns "chromosome,strand,full_position"; 
                __repr__ returns "Insertion_position(chromosome,strand,full_position)"
     - there are also _make_immutable and _make_mutable methods half-defined right now, for rudimentary and UNSAFE
                            hashability - those may be removed in the future (or implemented fully if there is need)
    """

    # NOTE: originally I used biopython SeqFeature objects for position_before and position_after (specifically SeqFeature.ExactPosition(min_val) for exactly positions and SeqFeature.WithinPosition(min_val, min_val-max_val) for ambiguous ones), but then I realized I'm not using that and it's over-complicated and hard to make immutable and may not be what I need anyway even if I do need immutable positions, so I switched to just integers. The function to generate those was called make_position_range, and I removed it on 2012-04-23, if I ever want to look it up again later.

    def __init__(self, chromosome, strand, full_position=None, position_before=None, position_after=None):
        """ Initialize all values - chromosome/strand are just copied from arguments; positions are more complicated. 
        
        You must provide either full_position, OR one or both of position_before/position_after. 
        The two position_* arguments must be castable to ints, or None.
        The full_position argument must be a string of the form '100-200', '?-200' or '100-?', such as would be generated 
         by self.full_position() - self.position_before and _after are set based on the two parts of the string.
        Self.min_/max_position are calculated based on self.position_before/_after - both, or whichever one isn't None.
        """
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
        return ','.join([self.chromosome, self.strand, self.full_position])

    def __repr__(self):
        return "Insertion_position('%s','%s','%s')"%(self.chromosome, self.strand, self.full_position)

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
    #  and since I also need a sane comparison operator, I can't use the default object id-based hashing, 
    #  so I have to make the objects immutable for them to be hashable. 
    # More info on how/why that is: http://docs.python.org/reference/datamodel.html#object.__hash__ 
    # Some links on implementation: http://stackoverflow.com/questions/9997176/immutable-dictionary-only-use-as-a-key-for-another-dictionary, http://stackoverflow.com/questions/1151658/python-hashable-dicts, http://stackoverflow.com/questions/4996815/ways-to-make-a-class-immutable-in-python, http://stackoverflow.com/questions/4828080/how-to-make-an-immutable-object-in-python
    # This implementation is not perfect, but it should do.

    # MAYBE-TODO implement!  DO I actually need hashability at all?  Right now the only time I put positions in sets is in mutant_join_datasets.py, which will be removed anyway... Will I want to do it again somewhere else? WAIT AND SEE.
    # MAYBE-TODO if I implement this properly, add to class docstring, and update method docstrings!
    # MAYBE-TODO if I implement this properly, remove the leading _ from the two _make* method names, since in that case they won't really be private any more - right now it's there just to remind me I shouldn't be using them.
    # MAYBE-TODO if I implement this properly, add to unit-tests!
    # MAYBE-TODO if I implement this properly, could add a hashable switch to __init__? And include it in __repr__.

    def _hash(self):
        """ Private hash-method based on _make_key - __hash__ will use it if you run self._make_immutable()."""
        return hash(self._make_key())

    def __hash__(self):
        """ If self.hashable is True, use private _hash method, otherwise use id-based object.__hash__. """
        if hasattr(self,'hashable') and self.hashable:  return self._hash()
        else:                                           return object.__hash__(self)

    def _make_immutable(self):
        """ Reversibly make object immutable and hashable. NOT FULLY IMPLEMENTED."""
        # MAYBE-TODO implement removing mutability methods!
        # set object to hashable - it's immutable, so that's all right
        self.hashable = True

    def _make_mutable(self):
        """ Reversibly make object mutable and non-hashable. NOT FULLY IMPLEMENTED."""
        # MAYBE-TODO implement re-adding mutability methods!
        # set object to unhashable, since it's mutable
        self.hashable = False


def get_insertion_pos_from_HTSeq_read_pos(HTSeq_pos, cassette_end, reads_are_reverse=False):
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
    return Insertion_position(chrom, strand, position_before=pos_before, position_after=pos_after)


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
    assert insertion_pos.strand in ['+','-'], "Strand should be %s, and is %s!"%(' or '.join(SEQ_STRANDS),strand)
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
            if gene.strand==1:      orientation = 'sense' if insertion_pos.strand=='+' else 'antisense'
            elif gene.strand==-1:   orientation = 'sense' if insertion_pos.strand=='-' else 'antisense'
            else:                   orientation = '?'
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


### Main classes describing the mutants and mutant sets

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
          want to apply the method to.  Doing this wrong will raise an Exception.  
         Some methods (like get_main_sequence) have a well-defined behavior when called on a multi-dataset mutant without 
          specifying a dataset - they'll give the result for all datasets added together.
         Some methods may not have multi-dataset functionality implemented, if I didn't think it would be useful.
    """

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
            self.by_dataset = defaultdict(lambda: Insertional_mutant(readcount_related_only=True))

    # MAYBE-TODO give each mutant some kind of unique ID at some point in the process?  Or is genomic location sufficient?  If we end up using per-mutant barcodes (in addition to the flanking sequences), we could use that, probably, or that plus genomic location.

    def _set_readcount_related_data_to_zero(self):
        """ Set all readcount-related data to 0/empty."""
        self.total_read_count      = 0
        self.perfect_read_count    = 0
        self.unique_sequence_count = 0
        self.sequences_and_counts  = defaultdict(lambda: 0)

    def _check_consistent_multi_dataset_args(self, dataset_name, function_name='<some_func>', multi_noname_allowed=False):
        """ Make sure dataset_name is consistent with single-/multi-dataset status of mutant; raise Exception if not.

        Raise Exception if dataset_name is None for multi-dataset mutant (i.e. the caller is trying to do some operation 
          on a multi-dataset mutant without specifying which dataset to do it on (there are a few cases in which this
           does make sense - in that case multi_noname_allowed=True should be passed), 
          or if dataset_name is not None for single-dataset mutant (i.e. the caller is trying to do some operation 
           on a specific dataset when the mutant is single-dataset and has no named datasets at all).
        The exception text will use function_name as the name of the calling function, if given.
        """
        # MAYBE-TODO probably possible to do this as a decorator and without passing function_name explicitly, but I'm not sure how...  See http://stackoverflow.com/questions/5063607/is-there-a-self-flag-can-reference-python-function-inside-itself for how to get the function name, but I don't actually know how to do the rest with a decorator (specifically how to get the dataset_name value from the decorated function so I can do things with it!)
        if self.multi_dataset and (dataset_name is None) and not multi_noname_allowed: 
            raise Exception("This mutant is in multi-dataset form! Provide dataset_name "
                            +"to specify which dataset to apply %s() to."%function_name)
        if not self.multi_dataset and (dataset_name is not None): 
            raise Exception("You're trying to apply %s() to dataset %s, but this mutant "%(function_name, dataset_name)
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

    def merge_mutant(self, other, dont_merge_positions=False, check_gene_data=True):
        """ Merge other mutant into this mutant: merge counts, sequences, optionally position; set other's counts to 0."""

        if self.multi_dataset:  
            raise Exception("Merging multi-dataset mutants not implemented!")
            # MAYBE-TODO implement merging for multi-dataset mutants?  Seems like unnecessary complication for now.

        # make sure the two mutants don't have conflicting gene data, if required (note: NOT checking position)
        #  (this needs to be done first, BEFORE we make changes to the mutants!)
        if check_gene_data:
            try:
                assert len(set([self.gene, other.gene]) - set([SPECIAL_GENE_CODES.not_determined])) <= 1
                assert len(set([self.orientation, other.orientation]) - set(['?'])) <= 1
                assert len(set([self.gene_feature, other.gene_feature]) - set(['?'])) <= 1
            except AssertionError:
                raise Exception("Can't merge the two mutants: the gene/orientation/feature data differs!")
        
        # merge read counts
        self.total_read_count += other.total_read_count
        # note: not incrementing perfect read counts, because "perfect" read counts from the wrong position aren't perfect!
        # merge sequences
        # LATER-TODO may want to keep more data about sequences! Like exact position and number of mutations - may want to store a list of HTSeq.alignment objects instead of just sequences+counts, really
        for (seq,count) in other.sequences_and_counts.iteritems():
            self.sequences_and_counts[seq] += count
        self.unique_sequence_count = len(self.sequences_and_counts)
        other._set_readcount_related_data_to_zero()
        # merge positions - MAYBE-TODO implement position merging?  Are we ever going to need it, really?
        if not dont_merge_positions:  
            raise Exception("Mutant merging with merging positions NOT IMPLEMENTED!")

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
        # Exception if dataset_name given for single-dataset mutant, but not giving it for multi-dataset is allowed
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

        ### now actually find the Nth main sequence in the chosen seq:count dictionary
        sequences_by_count = sorted([(count,seq) for (seq,count) in seqs_to_counts.iteritems()], reverse=True)
        # try returning the Nth sequence and count; return nothing if there are under N sequences.
        try:                return (sequences_by_count[N-1][1], sequences_by_count[N-1][0])
        except IndexError:  return ('',0)

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
                raise Exception("Can't run convert_to_multi_dataset - mutant is already multi-dataset!")
        # generate new by-dataset dictionary of all the readcount-related data separately for each dataset
        self.by_dataset = defaultdict(lambda: Insertional_mutant(readcount_related_only=True))
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
                                    force=False, check_constant_data=False):
        """ Copy all readcount-related data from other_mutant to self.by_dataset dictionary[other_mutant_dataset_name].

        If self isn't a multi-dataset mutant, raise an exception.
        If check_constant_data is True, check that the position/gene data of self and other_mutant matches.
        If self already has a other_mutant_dataset_name dataset, raise Exception, unless force=True, then overwrite it.
        """
        if not self.multi_dataset:
            raise Exception("This mutant hasn't been converted to multi-dataset form, "
                            +"can't run add_other_mutant_as_dataset! Run convert_to_multi_dataset first.")
        if other_mutant_dataset_name in self.by_dataset and not force:
            raise Exception("This mutant already has a %s dataset! Can't overwrite it with. "%other_mutant_dataset_name
                            +"new one.  Choose a different name for new dataset, or use force=True argument.")

        # if desired, check that the position/gene data matches 
        #  (probably should be using ifs rather than asserts, but I think since they're wrapped in a try/except it's fine)
        if check_constant_data:
            try:
                # MAYBE-TODO should this be an inexact comparison? Right now 100-101 is NOT equal to ?-101 or to 100-102.
                assert self.position == other_mutant.position
                assert len(set([self.gene, other_mutant.gene]) - set([SPECIAL_GENE_CODES.not_determined])) == 1
                assert len(set([self.orientation, other_mutant.orientation]) - set(['?'])) == 1
                assert len(set([self.gene_feature, other_mutant.gene_feature]) - set(['?'])) == 1
            except AssertionError:
                raise Exception("Can't add mutant2 as dataset to mutant1: the mutant position/gene data differs!")

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
            raise Exception("This mutant is not multi-dataset, can't run give_single_dataset_mutant!")
        if single_dataset_name not in self.by_dataset.keys() and not force:
            raise Exception("This mutant doesn't have a %s dataset! "%single_dataset_name
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
            raise Exception("This mutant is not multi-dataset, can't run give_all_single_dataset_mutants!")
        mutant_dictionary = defaultdict(lambda: Insertional_mutant(multi_dataset=False,readcount_related_only=False))
        for dataset_name in self.by_dataset.iterkeys():
            mutant_dictionary[dataset_name] = self.give_single_dataset_mutant(dataset_name)
        return mutant_dictionary

    # MAYBE-TODO should there be a copy_mutant function to make a deepcopy?


class Insertional_mutant_pool_dataset():
    """ A dataset of insertional mutants - contains some total counts and a set of Insertional_mutant objects.
    
    Attributes:
     - cassette_end - specifies which end of the insertion cassette the reads are on, and 
     - reads_are_reverse - True if the reads are in reverse orientation to the cassette, False otherwise
     - mutants_by_position - a (chrom,strand,pos):Insertional_mutant dictionary
     - discarded_read_count - number of reads discarded in preprocessing before alignment (not counted in total_read_count)
     - ignored_region_read_counts - region_name:read_count dictionary (not counted in total_read_count)
     - total_read_count, aligned_read_count, unaligned_read_count, perfect_read_count - various read counts, obvious
     - strand_read_counts, specific_region_read_counts - name:count dictionaries for strands and specific regions to track
     - mutants_in_genes, mutants_not_in_genes, mutants_undetermined - counts of mutants in genes, not in genes, unknown
     - mutant_counts_by_orientation, mutant_count_by_feature - name:count dictionaries for mutant gene location properties
    """

    # Implement new functionality for datasets:
    # - TODO joining multiple datasets and printing the output like in mutant_join_datasets.py (which I will rewrite)
    # - TODO subtracting datasets by presence (for contamination - remove ALL mutants from D1 that are present in D2, 
    #       or that have at least X reads in D2, or some such)
    # - MAYBE-TODO splitting joint datasets into separate ones?
    # - MAYBE-TODO adding and subtracting datasets (by readcount) - do we currently need that?
    # - TODO calculating growth rates!
    #
    # TODO once the multi-dataset version of Insertional_mutant_pool_dataset is finished and tested and works, implement dataset-joining using that instead of mutant_join_datasets.py!

    def __init__(self, cassette_end='?', reads_are_reverse='?'):
        """ Saves the arguments as properties of the dataset; initializes an empty mutant dict; sets all counts to 0. """
         # make sure the arguments are valid values
        if not cassette_end in SEQ_ENDS+['?']: 
            raise ValueError("The cassette_end variable must be one of %s or '?'!"%SEQ_ENDS)
        self.cassette_end = cassette_end
        if not reads_are_reverse in [True,False,'?']: 
            raise ValueError("The reads_are_reverse variable must be True, False, or '?'!")
        self.reads_are_reverse = reads_are_reverse
        # MAYBE-TODO should cassette_end and reads_are_reverse be specified for the whole dataset, or just for each set of data added, in add_alignment_reader_to_data? The only real issue with this would be that then I wouldn't be able to print this information in the summary - or I'd have to keep track of what the value was for each alignment reader added and print that in the summary if it's a single value, or 'varied' if it's different values. Might also want to keep track of how many alignment readers were involved, and print THAT in the summary!  Or even print each (infile_name, cassette_end, reads_are_reverse) tuple as a separate line in the header.
        # mutants_by_position is the main data structure here
        self.mutants_by_position = {}
        self.mutants_by_gene = defaultdict(lambda: [])
        # various total read/mutant counts to keep track of
        self.discarded_read_count = 'unknown'
        self.ignored_region_read_counts = defaultdict(lambda: 0)
        self.total_read_count, self.aligned_read_count, self.unaligned_read_count, self.perfect_read_count = 0,0,0,0
        self.strand_read_counts = defaultdict(lambda: 0)
        self.specific_region_read_counts = defaultdict(lambda: 0)
        self.mutants_in_genes, self.mutants_not_in_genes, self.mutants_undetermined = 0,0,0
        self.mutant_counts_by_orientation = defaultdict(lambda: 0)
        self.mutant_counts_by_feature = defaultdict(lambda: 0)
        self.mutant_merging_info = []
    
    def add_discarded_reads(self, N, reset_count=False):
        """ Add N to self.discarded_read_count (or set self.discarded_read_count to N if reset_count is True). """
        # LATER-TODO if preprocessing is changed to include more than one step at which reads can be discarded, make sure to keep track of that here and save each count separately! 
        if reset_count or self.discarded_read_count=='unknown':
            self.discarded_read_count = int(N)
        else:
            self.discarded_read_count += int(N)

    def add_alignment_reader_to_data(self, HTSeq_alignment_reader, uncollapse_read_counts=False, 
                                     treat_unknown_as_match=False, chromosomes_to_count=[], chromosomes_to_ignore=[]):
        """ Adds all alignments from the reader to the mutant data; currently based only on position, but that may change. 

        Input must be a list/generator/etc of HTSeq.Alignment objects (usually an HTSeq.SAM_Reader).
        Set uncollapse_read_counts to True if the original deepseq data was collapsed to unique sequences using
         fastx_uncollapser before alignment, to get the correct original read counts.
        Treat_unknown_as_match governs whether alignments with no detailed information are treated as perfect or not.
        Chromosomes_to_count is a list of chromosomes that should have aligned read counts kept and added to the header 
         summary (they're treated normally otherwise); reads that align to a chromosome in the chromosomes_to_ignore list 
         will be ignored in the data (but not the total counts contained in the header). 
        """

        if self.cassette_end == '?':
            raise Exception("Cannot add data from an alignment reader if cassette_end isn't specified! Please set the "
                +"cassette_end attribute of this Insertional_mutant_pool_dataset instance to one of %s first."%SEQ_ENDS)
        if self.reads_are_reverse == '?':
            raise Exception("Cannot add data from an alignment reader if reads_are_reverse isn't set! Please set the "
                +"reads_are_reverse attribute of this Insertional_mutant_pool_dataset instance to True/False first.")
        for aln in HTSeq_alignment_reader:
            if uncollapse_read_counts:      read_count = get_seq_count_from_collapsed_header(aln.read.name)
            else:                           read_count = 1
            self.total_read_count += read_count
            # if read is unaligned, add to unaligned count and skip to the next read
            if (not aln.aligned) or (aln.iv is None):
                self.unaligned_read_count += read_count
                continue
            # get the cassette insertion position (as an Insertion_position object)
            insertion_pos = get_insertion_pos_from_HTSeq_read_pos(aln.iv, self.cassette_end, self.reads_are_reverse)
            chrom = insertion_pos.chromosome
            strand = insertion_pos.strand
            pos = insertion_pos.min_position
            # if read is aligned to one of the chromosomes_to_ignore, add to the right count and skip to the next read
            if chrom in chromosomes_to_ignore:
                self.ignored_region_read_counts[chrom] += read_count
                continue
            # if read is aligned to anything else, add to aligned count, strand counts etc
            self.aligned_read_count += read_count
            self.strand_read_counts[strand] += read_count
            # MAYBE-TODO do I want info on how many reads were aligned to which strand for the chromosomes_to_count or even all chromosomes?  Maybe optionally...  And how many were perfect, and how many mutants there were... Might want to write a separate class or function just for this.  If so, should output it in a tabular format, with all the different data (reads, +, -, perfect, ...) printed tab-separated in one row.
            if chrom in chromosomes_to_count:
                self.specific_region_read_counts[chrom] += read_count
            # check if there's already a mutant at that position; if not, make a new one
            if (chrom,strand,pos) not in self.mutants_by_position.keys():
                self.mutants_by_position[(chrom,strand,pos)] = Insertional_mutant(insertion_pos)
            # add_read adds the read to the full data; also returns True if alignment was perfect, to add to perfect_count
            if self.mutants_by_position[(chrom,strand,pos)].add_read(aln, read_count, 
                                                                     treat_unknown_as_match=treat_unknown_as_match):
                self.perfect_read_count += read_count


    def merge_adjacent_mutants(self, merge_max_distance=1, merge_count_ratio=1000, dont_change_positions=False, verbosity_level=1):
        """ Merge adjacent mutants based on strand, distance, and count ratio; optionally merge positions; print counts.
        Merge mutants if they're distant by merge_max_distance or less and one has at least merge_count_ratio x fewer 
         reads than the other; merge their read counts, sequences etc, and remove the lower-read-count one from 
         self.mutants_by_position.  
        If dont_change_positions is True, also change the full position data of the remaining mutant to reflect the merge.
        """

        # MAYBE-TODO After this merging, a dictionary by position would no longer really make sense, would it?  Well, I guess it sort of would, really.  Should be fine for now.  But I could change it to just a list - the dictionary by positon isn't used that often, except in mutant_join_datasets.py, which needs rewriting anyway.

        if verbosity_level>1:
            print "Merging adjacent mutants: max distance %s, min count ratio %s..."%(merge_max_distance,merge_count_ratio)
        all_positions = sorted(self.mutants_by_position.keys())
        adjacent_opposite_strands_count = 0
        adjacent_same_strands_merged = 0
        adjacent_same_strands_left = 0
        # go over each pos1,pos2 pair (only once)
        for ((chr1,strand1,pos1), (chr2,strand2,pos2)) in combinations(all_positions,2):
            # if the two positions are on different chromosomes or aren't adjacent, skip
            if not (chr1==chr2 and abs(pos2-pos1)<=merge_max_distance):  
                continue
            # if the two positions are adjacent but on different strands, they can't be merged, count and skip
            if not strand1==strand2:
                adjacent_opposite_strands_count += 1
                if verbosity_level>1:
                    print "  LEAVING adjacent mutants: %s strand %s %s and strand %s %s."%(chr1,strand1,pos1,strand2,pos2)
                continue
            # make sure these positions still exist in self.mutants_by_position (haven't been merged)
            try:
                mutant1 = self.mutants_by_position[(chr1,strand1,pos1)]
            except KeyError:
                if verbosity_level>1:
                    print "Warning: attempting to merge a mutant that no longer exists %s with %s!"%((chr1,strand1,pos1),
                                                                                                     (chr2,strand2,pos2))
                continue
            try:
                mutant2 = self.mutants_by_position[(chr2,strand2,pos2)]
            except KeyError:
                if verbosity_level>1:
                    print "Warning: attempting to merge a mutant that no longer exists %s with %s!"%((chr2,strand2,pos2),
                                                                                                     (chr1,strand1,pos1))
                continue
            mutant1_readcount, mutant2_readcount = mutant1.total_read_count, mutant2.total_read_count
            readcount_ratio = max(mutant1_readcount/mutant2_readcount, mutant2_readcount/mutant1_readcount)
            # if the two mutants are adjacent, but the readcounts aren't different enough , count and skip
            if not readcount_ratio>=merge_count_ratio:
                adjacent_same_strands_left += 1
                if verbosity_level>1:
                    print "  LEAVING adjacent mutants: %s strand %s %s and strand %s %s,"%(chr1,strand1,pos1,strand2,pos2),
                    print "%s and %s reads."%(mutant1_readcount, mutant2_readcount)
                continue
            # if the two mutants are adjacent and with sufficiently different readcounts, count, 
            #  MERGE the lower-readcount mutant into the higher-readcount one, 
            #  and remove the lower one from the mutants_by_position dictionary
            adjacent_same_strands_merged += 1
            if verbosity_level>1:
                print " MERGING adjacent mutants: %s, strand %s %s and strand %s %s,"%(chr1,strand1,pos1,strand2,pos2),
                print "%s and %s reads."%(mutant1_readcount, mutant2_readcount)
            if mutant1_readcount > mutant2_readcount:
                mutant1.merge_mutant(mutant2, dont_merge_positions=dont_change_positions)
                del self.mutants_by_position[(chr2,strand2,pos2)]
            elif mutant2_readcount > mutant1_readcount:
                mutant2.merge_mutant(mutant1, dont_merge_positions=dont_change_positions)
                del self.mutants_by_position[(chr1,strand1,pos1)]
        self.mutant_merging_info.append("merged %s pairs of adjacent mutants "%adjacent_same_strands_merged 
                                        +"(max distance %s, min count ratio %s) "%(merge_max_distance, merge_count_ratio))
        self.mutant_merging_info.append("  (kept %s pairs on same strand "%adjacent_same_strands_left 
                                        +"and %s on opposite strands)"%adjacent_opposite_strands_count)
        # LATER-TODO add unit-tests or run-tests for this!


    def find_genes_for_mutants(self, genefile, detailed_features=False, known_bad_chromosomes=[], 
                               N_run_groups=3, verbosity_level=1):
        """ To each mutant in the dataset, add the gene it's in (look up gene positions for each mutant using genefile).

        If detailed_features is True, also look up whether the mutant is in an exon/intron/UTR (NOT IMPLEMENTED); 
        Read the file in N_run_groups passes to avoid using up too much memory/CPU.
        """ 

        # group all the mutants by chromosome, so that I can go over each chromosome in genefile separately
        #   instead of reading in all the data at once (which uses a lot of memory)
        mutants_by_chromosome = defaultdict(lambda: set())
        for (chrom,_,_),mutant in self.mutants_by_position.iteritems():
            mutants_by_chromosome[chrom].add(mutant)

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
        for chromosome_set in chromosome_sets:
            genefile_parsing_limits = {'gff_id': list(chromosome_set)}
            if not detailed_features: 
                genefile_parsing_limits['gff_type'] = ['gene']
            with open(genefile) as GENEFILE:
                for chromosome_record in GFF.parse(GENEFILE, limit_info=genefile_parsing_limits):
                    if verbosity_level>1:   print "    parsing %s for mutant gene locations..."%chromosome_record.id
                    for mutant in mutants_by_chromosome[chromosome_record.id]:
                        gene_ID, orientation, feature = find_gene_by_pos(mutant.position, chromosome_record, 
                                                                         detailed_features, quiet=(verbosity_level==0))
                        mutant.gene, mutant.orientation, mutant.gene_feature = gene_ID, orientation, feature
                        if gene_ID==SPECIAL_GENE_CODES.not_found:   self.mutants_not_in_genes += 1
                        else:                                       self.mutants_in_genes += 1
                        if orientation not in ['?','-']:            self.mutant_counts_by_orientation[orientation] += 1
                        if feature not in ['?','-']:                self.mutant_counts_by_feature[feature] += 1

        # for mutants in chromosomes that weren't listed in the genefile, use special values
        for chromosome in set(mutants_by_chromosome.keys())-set(all_reference_chromosomes):
            if not chromosome in known_bad_chromosomes:
                print 'Warning: chromosome "%s" not found in genefile data!'%(chromosome)
            for mutant in mutants_by_chromosome[chromosome]:
                mutant.gene,mutant.orientation,mutant.gene_feature = SPECIAL_GENE_CODES.chromosome_not_in_reference,'-','-'
                self.mutants_undetermined += 1

    def make_by_gene_mutant_dict(self):
        """ Fill the self.mutants_by_gene dictionary (gene_name:mutant_list) based on self.mutants_by_position; 
        real gene IDs only (ignore SPECIAL_GENE_CODES).
        """
        for mutant in self.mutants_by_position.itervalues():
            if mutant.gene not in SPECIAL_GENE_CODES.all_codes:
                self.mutants_by_gene[mutant.gene].append(mutant)
        # LATER-TODO add unit-test

    def gene_dict_by_mutant_number(self):
        """ Return a mutant_count:gene_ID_set dictionary of genes with 1/2/etc mutants. """
        self.make_by_gene_mutant_dict()
        gene_by_mutantN = defaultdict(lambda: set())
        for (gene,mutants) in self.mutants_by_gene.iteritems():
            gene_by_mutantN[len(mutants)].add(gene)
        # check that the numbers add up to the total number of genes
        all_genes = set([m.gene for m in self.mutants_by_position.itervalues()]) - set(SPECIAL_GENE_CODES.all_codes)
        assert sum([len(geneset) for geneset in gene_by_mutantN.itervalues()]) == len(all_genes)
        return gene_by_mutantN
        #LATER-TODO add unit-test

    def add_detailed_gene_info(self, gene_info_file):
        """ _____ """
        print "    *** WARNING: getting detailed gene info (-G option) NOT IMPLEMENTED! Skipping. ***"
        # TODO implement!  Already done in mutant_join_datasets.py...  Get gene name/description/stuff from gene_info_file and add to mutants (how, exactly?  Do we want to add those fields to each mutant separately, or make a new dictionary of genes that can be looked up by geneID to get the detailed info, or just make gene objects that contain all this stuff and use those as mutant.gene instead of just geneID strings?
        # Once it's implemented, make sure it gets printed somewhere - probably only in the line-per-gene output
        # LATER-TODO add this to the gene-info run-test case!

    def find_most_common_mutant(self):
        """ Return the Insertional_mutant object from self.mutants_by_position with the most total reads."""
        current_readcount = 0
        current_mutant = None
        for mutant in self.mutants_by_position.itervalues():
            if mutant.total_read_count > current_readcount:
                current_readcount = mutant.total_read_count
                current_mutant = mutant
        return current_mutant

    def print_summary(self, OUTPUT=sys.stdout, N_genes_to_print=5, line_prefix='    ', header_prefix=' * '):
        """ Print basic read and mutant counts (prints to stdout by default, can also pass an open file object)."""
        # TODO should probably rewrite this to start with a real total read number (discarded+total)
        OUTPUT.write(line_prefix+"Reads discarded in preprocessing: %s\n"%(self.discarded_read_count))
        OUTPUT.write(header_prefix+"Total reads processed: %s\n"%(self.total_read_count))
        OUTPUT.write(line_prefix+"Unaligned reads: %s\n"%(self.unaligned_read_count))
        for (region,count) in self.ignored_region_read_counts.iteritems():
            OUTPUT.write(line_prefix+"Discarded reads aligned to %s: %s\n"%(region, count))
        OUTPUT.write(line_prefix+"Aligned reads (non-discarded): %s\n"%(self.aligned_read_count))
        OUTPUT.write(line_prefix+"Perfectly aligned reads (no mismatches): %s\n"%(self.perfect_read_count))
        for (strand,count) in self.strand_read_counts.iteritems():
            OUTPUT.write(line_prefix+"Reads with insertion direction matching chromosome %s strand: %s\n"%(strand, count))
        for (region,count) in sorted(self.specific_region_read_counts.iteritems()):
            OUTPUT.write(line_prefix+"Reads aligned to %s: %s\n"%(region, count))
        # TODO add percentages of total (or aligned) reads to all of these numbers in addition to raw counts!
        # MAYBE-TODO keep track of the count of separate mutants in each category, as well as total read counts?
        OUTPUT.write(header_prefix+"Distinct mutants (read groups) by cassette insertion position: %s\n"%\
                                                                                     (len(self.mutants_by_position)))
        OUTPUT.write(line_prefix+"(read is at %s end of cassette, in %s direction to cassette)\n"%(self.cassette_end, 
                                                {'?': '?', True: 'reverse', False: 'forward'}[self.reads_are_reverse]))
        for line in self.mutant_merging_info:
            OUTPUT.write(line_prefix+line+'\n')
        g = self.find_most_common_mutant()
        OUTPUT.write(line_prefix+"Most common mutant: %s, position %s, %s strand:"%(g.position.chromosome,
                                                                           g.position.full_position, g.position.strand))
        OUTPUT.write(" %s reads (%.0f%% of aligned)\n"%(g.total_read_count, 
                                                        100.0*g.total_read_count/self.aligned_read_count))
        # MAYBE-TODO may also be a good idea to keep track of the most common SEQUENCE, not just mutant...
        # print the gene annotation info, but only if there is any
        if self.mutants_in_genes + self.mutants_not_in_genes + self.mutants_undetermined:
            OUTPUT.write(line_prefix+"Mutant cassettes in unknown chromosomes: %s\n"%(self.mutants_undetermined))
            OUTPUT.write(line_prefix+"Mutant cassettes not inside genes: %s\n"%(self.mutants_not_in_genes))
            OUTPUT.write(header_prefix+"Mutant cassettes inside genes: %s\n"%(self.mutants_in_genes))
            for (orientation,count) in sorted(self.mutant_counts_by_orientation.items(),reverse=True):
                OUTPUT.write(line_prefix+"Mutant cassettes in %s orientation to gene: %s\n"%(orientation,count))
            # custom order for features to make it easier to read: CDS, intron, UTRs, everything else alphabetically after
            proper_feature_order = defaultdict(lambda: 3, {'CDS':0, 'intron':1, 'five_prime_UTR':2, 'three_prime_UTR':2})
            for (feature,count) in sorted(self.mutant_counts_by_feature.items(), 
                                          key=lambda (f,n): (proper_feature_order[f],f)):
                OUTPUT.write(line_prefix+"Mutant cassettes in gene feature %s: %s\n"%(feature,count))
            all_genes = set([m.gene for m in self.mutants_by_position.itervalues()]) - set(SPECIAL_GENE_CODES.all_codes)
            OUTPUT.write(header_prefix+"Genes containing a mutant: %s\n"%(len(all_genes)))
            genes_by_mutantN = self.gene_dict_by_mutant_number()
            for (mutantN, geneset) in genes_by_mutantN.iteritems():
                if N_genes_to_print>0:
                    genelist_to_print = ', '.join(sorted(list(geneset))[:N_genes_to_print])
                    if len(geneset)<=N_genes_to_print:  genelist_string = ' (%s)'%genelist_to_print
                    else:                               genelist_string = ' (%s, ...)'%genelist_to_print
                else:                                   genelist_string = ''
                OUTPUT.write(line_prefix+"Genes with %s mutants: %s%s\n"%(mutantN, len(geneset), genelist_string))
            # TODO Add some measure of mutations, like how many mutants have <50% perfect reads (or something - the number should probably be a command-line option).  Maybe how many mutants have <20%, 20-80%, and >80% perfect reads (or 10 and 90, or make that a variable...)


    def print_data(self, OUTPUT=None, sort_data_by=None, N_sequences=2, header_line=True, header_prefix="# "):
        """ Print full data, one line per mutant: position data, gene info, read counts, optionally sequences.
        (see the file header line for exactly what all the output fields are).

        If N_sequences<1, only prints position (chromosome and position), total and perfect read count, 
          and the number of sequence variants.  If N_sequences==1, also prints the most common sequence and count; 
          if N_sequences>1, adds the second most common sequence and count, and so on.
        Output is tab-separated, with headers starting with "# ".  Prints to an open file object (stdout by default).
        """
        # MAYBE-TODO should printing the gene info be optional?  Maybe... Would save space when there isn't any meaningful gene info to print anyway.
        if OUTPUT is None:
            OUTPUT = sys.stdout

        if header_line:
            headers = ['chromosome','strand','min_position','full_position', 'gene','orientation','feature',
                       'total_reads','perfect_reads', 'N_sequence_variants']
            for N in range(1,N_sequences+1):
                headers.extend(['read_sequence_%s'%N,'seq_%s_count'%N])
            OUTPUT.write(header_prefix + '\t'.join(headers) + "\n")

        if sort_data_by=='position':
            data = sorted(self.mutants_by_position.values(), key = lambda x: x.position)
            # x.position here is an Insertion_position object and has a sensible cmp function
        elif sort_data_by=='read_count':
            data = sorted(self.mutants_by_position.values(), 
                          key = lambda x: (x.total_read_count, x.perfect_read_count, x.position), reverse=True)
        else:
            data = self.mutants_by_position.itervalues()

        for mutant in data:
            mutant_data = [mutant.position.chromosome, mutant.position.strand, mutant.position.min_position, 
                           mutant.position.full_position, mutant.gene, mutant.orientation, mutant.gene_feature, 
                           mutant.total_read_count, mutant.perfect_read_count, mutant.unique_sequence_count]
            OUTPUT.write('\t'.join([str(x) for x in mutant_data]))
            for N in range(1,N_sequences+1):
                OUTPUT.write('\t%s\t%s'%mutant.get_main_sequence(N))
                # MAYBE-TODO also give the length and number of mutations for each sequence? Optionally?  Length is easy, but do I even keep track of mutation number?  I probably should...
            OUTPUT.write('\n')

    def read_from_file(self, infile, assume_new_sequences=False):
        """ Read data from a file generated by self.print_data, and add to self.mutants_by_position. Ignores some things. 
        
        Populates most of the dataset total read/mutant count values correctly, but ignores unaligned and discarded reads
         and anything involving special regions.
        Ignores single sequence/count fields; the total number of sequence variants is unreliable if you add to preexisting
        data (if the data originally listed 1 unique sequence, and new data adds another 2 sequences, is the total 2 or 3 
         unique sequences?  If assume_new_sequences is True, the total is old+new; if it's False, it's max(old,new)). 
        """
        for line in open(infile):
            # LATER-TODO get unaligned/discarded/etc read count from summary, so we can keep track of full counts!
            # ignore comment and header lines, parse other tab-separated lines into values
            if line.startswith('#'):    continue
            if line.startswith('chromosome\tstrand\tmin_position\t'):    continue
            fields = line.split('\t')
            chromosome = fields[0]
            strand = fields[1]
            min_pos = int(fields[2])
            full_pos = fields[3]
            gene, orientation, gene_feature = fields[4:7]
            total_reads,perfect_reads,sequence_variants = [int(x) for x in fields[7:10]]
            # generate a new mutant if necessary; add the counts and gene info to the mutant
            basic_position_data = (chromosome,strand,min_pos)
            if basic_position_data not in self.mutants_by_position.keys():
                full_position_data = Insertion_position(chromosome, strand, full_position=full_pos)
                self.mutants_by_position[basic_position_data] = Insertional_mutant(full_position_data)
            self.mutants_by_position[basic_position_data].add_counts(total_reads, perfect_reads, 
                                                                           sequence_variants, assume_new_sequences)
            # MAYBE-TODO I'm overwriting old gene data here, even if the new data is '-'!  Write a better update function.
            self.mutants_by_position[basic_position_data].gene = gene
            self.mutants_by_position[basic_position_data].orientation = orientation
            self.mutants_by_position[basic_position_data].gene_feature = gene_feature
            # get however many specific sequences/counts are listed (this is variable)
            sequence_fields = fields[10::2]
            count_fields = fields[11::2]
            for seq, count in zip(sequence_fields, count_fields):
                if int(count)>0:
                    assert seq!=''
                    self.mutants_by_position[basic_position_data].sequences_and_counts[seq] += int(count)
            # add to dataset total read/mutant counts
            # MAYBE-TODO Might just want to write a function that recalculates all of the total counts below, to be ran at the end of read_from_file and add_alignment_reader_to_data and I guess find_genes_for_mutants, instead of doing it this way
            self.total_read_count += total_reads
            self.aligned_read_count += total_reads
            self.perfect_read_count += perfect_reads
            self.strand_read_counts[strand] += total_reads
            if gene==SPECIAL_GENE_CODES.not_found:        self.mutants_not_in_genes += 1
            elif gene in SPECIAL_GENE_CODES.all_codes:    self.mutants_undetermined += 1  # the two codes beside not_found
            else:                                         self.mutants_in_genes += 1
            if orientation not in ['?','-']:              self.mutant_counts_by_orientation[orientation] += 1
            if gene_feature not in ['?','-']:             self.mutant_counts_by_feature[gene_feature] += 1
    

# LATER-TODO at some point I'll probably want to generalize both of those classes to not necessarily be grouped by position...


######### Test code #########

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
        assert all([str(pos) == 'chr,+,100-?' for pos in all_before_positions])
        assert all([repr(pos) == "Insertion_position('chr','+','100-?')" for pos in all_before_positions])
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
        self.assertRaises(Exception, mutant.add_read, perfect_aln, read_count=3, dataset_name='d1')
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
        self.assertRaises(Exception, mutant.add_read, perfect_aln, read_count=3)

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
        mutant1.merge_mutant(mutant2, dont_merge_positions=True)
        assert (mutant1.total_read_count, mutant1.perfect_read_count) == (3, 2)
        assert (mutant2.total_read_count, mutant2.perfect_read_count) == (0, 0)
        assert (mutant1.sequences_and_counts, mutant2.sequences_and_counts) == ({'AAA':2, 'AA':1}, {})
        # if we define the gene for one of the mutants but not the other, still works
        #  (note that starting here mutant2 has 0 reads, but we can still merge it all we want just to check for errors)
        mutant2.gene = 'X'
        mutant1.merge_mutant(mutant2, check_gene_data=True, dont_merge_positions=True)
        # if we define the gene for both mutants and it's the same gene, also works
        mutant1.gene = 'X'
        mutant1.merge_mutant(mutant2, check_gene_data=True, dont_merge_positions=True)
        # if we define the gene for both mutants and it's a DIFFERENT gene, we get an exception unless we don't check
        mutant2.gene = 'Y'
        self.assertRaises(Exception, mutant1.merge_mutant, mutant2, dont_merge_positions=True, check_gene_data=True)
        mutant1.merge_mutant(mutant2, check_gene_data=False, dont_merge_positions=True)
        # dont_merge_positions=False option currently NOT IMPLEMENTED
        self.assertRaises(Exception, mutant1.merge_mutant, mutant2, dont_merge_positions=False, check_gene_data=False)
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
        self.assertRaises(Exception, multi_mutant.add_counts, 1,1,1)
        self.assertRaises(Exception, mutant.add_counts, 1,1,1, dataset_name='d1')

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
        self.assertRaises(Exception, multi_mutant.add_sequence_and_counts, 'GGG',1)
        self.assertRaises(Exception, mutant.add_sequence_and_counts, 'GGG',1, dataset_name='d1')

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
        # converting an already multi-dataset mutant to multi-dataset: Exception if not ignore_if_already_multi, 
        #  otherwise works and doesn't change anything - all the same asserts
        self.assertRaises(Exception, mutant1.convert_to_multi_dataset, current_dataset_name='d')
        self.assertRaises(Exception, mutant2.convert_to_multi_dataset)
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
        self.assertRaises(Exception, mutant1.add_other_mutant_as_dataset, mutant2, 'd2')
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
        # can't add overwrite an existing dataset name, unless force==True
        self.assertRaises(Exception, mutant1.add_other_mutant_as_dataset, mutant2, 'd2')
        mutant1.add_other_mutant_as_dataset(mutant2, 'd2', force=True)
        # if the two mutants have different positions, it should fail, unless check_constant_data=False
        mutant3 = Insertional_mutant(Insertion_position('chr','+',position_before=5))
        mutant4 = Insertional_mutant(Insertion_position('chr2','+',position_before=3))
        mutant5 = Insertional_mutant(Insertion_position('chr','-',position_before=3))
        self.assertRaises(Exception, mutant1.add_other_mutant_as_dataset, mutant3, 'd3', check_constant_data=True)
        self.assertRaises(Exception, mutant1.add_other_mutant_as_dataset, mutant4, 'd4', check_constant_data=True)
        self.assertRaises(Exception, mutant1.add_other_mutant_as_dataset, mutant5, 'd5', check_constant_data=True)
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
        self.assertRaises(Exception, mutant.give_single_dataset_mutant, 'd0', force=False)
        new_mutant_0 = mutant.give_single_dataset_mutant('d0', force=True)
        assert new_mutant_0.position == mutant.position
        assert new_mutant_0.gene == mutant.gene
        assert new_mutant_0.total_read_count == 0
        assert new_mutant_0.perfect_read_count == 0
        assert new_mutant_0.unique_sequence_count == 0
        assert new_mutant_0.sequences_and_counts == {}
        # trying to extract a dataset from a single-dataset mutant should fail
        mutant1 = Insertional_mutant(Insertion_position('chr','+',position_before=3), multi_dataset=False)
        self.assertRaises(Exception, mutant1.give_single_dataset_mutant, 'd2')

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
                assert data.cassette_end == cassette_end
                assert data.reads_are_reverse == reads_are_reverse
                assert data.mutants_by_position == {}
                assert data.total_read_count == data.aligned_read_count == data.perfect_read_count == 0
                assert data.unaligned_read_count == 0
                assert data.discarded_read_count == 'unknown'
                assert data.ignored_region_read_counts == data.strand_read_counts == data.specific_region_read_counts == {}
                assert data.mutants_in_genes == data.mutants_not_in_genes == data.mutants_undetermined == 0
                assert data.mutant_counts_by_orientation == data.mutant_counts_by_feature == {}
        for cassette_end in [True, False, None, 0, 1, 0.11, 23, 'asdfas', '', 'something', [2,1], {}]:
            self.assertRaises(ValueError, Insertional_mutant_pool_dataset, cassette_end, '?')
        for reads_are_reverse in ['forward', 'reverse', None, 0.11, 23, 'asdfas', '', 'something', [2,1], {}]:
            # note that this list doesn't include 0/1 because 0==False and 1==True, and that's what "in" tests 
            self.assertRaises(ValueError, Insertional_mutant_pool_dataset, '?', reads_are_reverse)

    # LATER-TODO add unit-test for add_discarded_reads, find_genes_for_mutants, find_most_common_mutant, 

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

    def test__read_from_file(self):
        ## 1. input file with no gene information but more variation in other features
        input_file = 'test_data/test_output__5prime.txt'
        data = Insertional_mutant_pool_dataset()
        data.read_from_file(input_file)
        assert data.total_read_count == data.aligned_read_count == 30
        assert data.unaligned_read_count == 0
        assert data.perfect_read_count == 22
        assert data.strand_read_counts == {'+':27, '-':3}
        assert len(data.mutants_by_position) == 17
        assert data.mutants_in_genes == data.mutants_not_in_genes == 0
        assert data.mutants_undetermined == 17
        assert data.mutant_counts_by_orientation == {}
        assert data.mutant_counts_by_feature == {}
        for mutant in data.mutants_by_position.itervalues():
            assert mutant.gene == SPECIAL_GENE_CODES.not_determined
            assert mutant.orientation == mutant.gene_feature == '?'
        # just spot-checking some of the outputs
        mutant = data.mutants_by_position[('reads_2_seqs_1','+',204)]
        assert mutant.position.chromosome == 'reads_2_seqs_1'
        assert mutant.position.min_position == 204
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 2
        assert mutant.perfect_read_count == 2 
        assert mutant.unique_sequence_count == 1
        mutant = data.mutants_by_position[('mutation_yes','+',204)]
        assert mutant.position.chromosome == 'mutation_yes'
        assert mutant.position.min_position == 204
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 6
        assert mutant.perfect_read_count == 0
        assert mutant.unique_sequence_count == 1
        mutant = data.mutants_by_position[('strandedness_-_normal_+_reverse','-',100)]
        assert mutant.position.chromosome == 'strandedness_-_normal_+_reverse'
        assert mutant.position.min_position == 100
        assert mutant.position.strand == '-'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1
        assert mutant.unique_sequence_count == 1
        ## 2. adding more data to a file that already has some...
        data.read_from_file(input_file, assume_new_sequences=False)
        mutant = data.mutants_by_position[('reads_2_seqs_1','+',204)]
        assert mutant.position.chromosome == 'reads_2_seqs_1'
        assert mutant.position.min_position == 204
        assert mutant.total_read_count == 4
        assert mutant.perfect_read_count == 4
        # how mutant.unique_sequence_count should act in this case depends on the value of assume_new_sequences
        assert mutant.unique_sequence_count == 1
        data.read_from_file(input_file, assume_new_sequences=True)
        assert mutant.unique_sequence_count == 2
        ## 3. input file with gene information
        input_file2 = 'test_data/test_output2__with-genes.txt'
        data2 = Insertional_mutant_pool_dataset()
        data2.read_from_file(input_file2)
        assert data2.total_read_count == data2.aligned_read_count == data2.perfect_read_count == 40
        assert data2.unaligned_read_count == 0
        assert data2.strand_read_counts == {'+':38, '-':2}
        assert len(data2.mutants_by_position) == 40
        assert data2.mutants_in_genes == 39
        assert data2.mutants_not_in_genes == 1
        assert data2.mutants_undetermined == 0
        assert data2.mutant_counts_by_orientation == {'sense':37, 'antisense':2}
        assert data2.mutant_counts_by_feature['CDS'] == 6
        assert data2.mutant_counts_by_feature['??'] == 2
        assert data2.mutant_counts_by_feature['CDS/three_prime_UTR'] == 1
        mutant = data2.mutants_by_position[('chromosome_A','+',20)]
        assert mutant.position.chromosome == 'chromosome_A'
        assert mutant.position.min_position == 20
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1 
        assert mutant.unique_sequence_count == 1
        assert mutant.gene == SPECIAL_GENE_CODES.not_found
        assert mutant.orientation == mutant.gene_feature == '-'
        mutant = data2.mutants_by_position[('chromosome_A','+',150)]
        assert mutant.position.chromosome == 'chromosome_A'
        assert mutant.position.min_position == 150
        assert mutant.position.strand == '+'
        assert mutant.total_read_count == 1
        assert mutant.perfect_read_count == 1 
        assert mutant.unique_sequence_count == 1
        assert mutant.gene == "test.geneA0_proper_plus"
        assert mutant.orientation == 'sense'
        assert mutant.gene_feature == 'five_prime_UTR'


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    # if I wanted more control I could do this instead:
    #import os
    #unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromName(os.path.splitext(sys.argv[0])[0]))
    #   (autodetection of all tests - see http://docs.python.org/library/unittest.html#unittest.TestLoader)
    # there's probably also some way to easily get all tests from the current file without passing the name, but I haven't found it yet...
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main(argv=[sys.argv[0]])

