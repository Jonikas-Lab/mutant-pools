#!/usr/bin/env python
""" Look at a gff file, print some overview info. 
 - print the basic structure of the records in the file (gene->mRNA->UTRs/exons or such)
 - print the gff IDs, sources, types etc and their counts
 - print the lengths (approximate) of each sequence, the number of genes it contains,
    and what percentage of the length is covered by genes
 - ...

WORK IN PROGRESS.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: gff_check_issues_coverage.py [options] gff_infile """

# basic libraries
import sys, time
import unittest
import itertools
from collections import defaultdict
import pprint
# other libraries
from BCBio import GFF
# my modules
from general_utilities import count_list_values, split_into_N_sets_by_counts


def check_for_overlapping_genes(sequence_record):
    """ Given a GFF sequence record (from BCBio.GFF parser), return list of tuples of IDs of overlapping genes.  
    Not perfect: may miss some overlaps if some genes are completely contained within other genes (raises warning), 
     but if it returns no overlaps, there definitely aren't any."""
    overlapping_gene_pairs = []
    all_gene_positions = []
    for gene in sequence_record.features:
        # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3), 
        #   so add 1 to start and keep end as is to convert to 1-based-end-inclusive
        all_gene_positions.append((gene.location.start.position+1, gene.location.end.position, gene.id))
    all_gene_positions.sort()
    for gene1_data,gene2_data in itertools.izip(all_gene_positions,all_gene_positions[1:]):
        (gene1_start,gene1_end,gene1_name), (gene2_start,gene2_end,gene2_name) = gene1_data, gene2_data
        if gene1_end>=gene2_start:
            overlapping_gene_pairs.append((gene1_name,gene2_name))
        # check for "gene1 contains gene2", print a warning, since it can make other things not work right
        if gene1_end>=gene2_end:
            print("WARNING: gene %s is completely inside gene %s! "%(gene1_name, gene2_name)
                      +"Various gene-position-related results may be inaccurate.")
    return overlapping_gene_pairs
    # MAYBE-TODO rewrite it so it actually detects ALL overlaps?  Right now if gene A contains nonoverlapping genes B and C, it'll sort them as (A,B,C) since A starts first, so it'll detect the (A,B) overlap, but it won't detect the (A,C) overlap because it doesn't CHECK (A,C), only (A,B) and (B,C).  This could be fixed either by just brute-force checking all gene pairs (and then using DNA_basic_utilities.position_test_overlap), or by writing something prettier.  In any case, not a priority, since generally genes DON'T OVERLAP AT ALL.


def check_for_overlapping_features(mRNA, gene_name):
    """ Given a GFF gene record (from BCBio.GFF parser), return list of tuples of IDs of overlapping features.  
    Not perfect: may miss some overlaps if some features are completely contained within other features, 
     but if it returns no overlaps, there definitely aren't any."""
    overlap_found = False
    all_feature_positions = []
    for feature in mRNA.sub_features:
        # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3), 
        #   so add 1 to start and keep end as is to convert to 1-based-end-inclusive
        all_feature_positions.append((feature.location.start.position+1, feature.location.end.position, feature.id))
    all_feature_positions.sort()
    for feature1_data,feature2_data in itertools.izip(all_feature_positions,all_feature_positions[1:]):
        (_,feature1_end,feature1_name), (feature2_start,feature2_end,feature2_name) = feature1_data, feature2_data
        if feature1_end>=feature2_start:
            print "WARNING: feature %s overlaps feature %s in gene %s! "%(feature1_name, feature2_name, gene_name)
            overlap_found = True
        # check for "feature1 contains feature2", print a warning, since it can make other things not work right
        if feature1_end>=feature2_end:
            print("WARNING: feature %s is completely inside feature %s in gene %s!"%(feature1_name,feature2_name,gene_name)
                      +" There may be additional undetected overlaps downstream.")
    return overlap_found


def find_min_gene_distance(sequence_record, starting_values=None):
    """ Given a GFF record from BCBio.GFF parser, return the lowest distance between adjacent genes, and the two gene IDs. 
    May not be accurate if some genes are completely contained inside other genes. """
    min_distance = sequence_record.seq._length if starting_values is None else starting_values[0]
    min_gene1 = 'none' if starting_values is None else starting_values[1]
    min_gene2 = 'none' if starting_values is None else starting_values[2]
    all_gene_positions = []
    for gene in sequence_record.features:
        # BCBio uses 0-based and end-exclusive positions (first-third base is bases 0,1,2, i.e range 0-3), 
        #   so add 1 to start and keep end as is to convert to 1-based-end-inclusive
        all_gene_positions.append((gene.location.start.position, gene.location.end.position-1, gene.id))
    all_gene_positions.sort()
    for (_,gene1_end,gene1_name), (gene2_start,_,gene2_name) in itertools.izip(all_gene_positions,all_gene_positions[1:]):
        # subtract 1 from distance, so if gene1 is 1-4 and gene2 is 5-9 the distance is 0
        gene_distance = gene2_start - gene1_end - 1
        if gene_distance < min_distance:
            min_distance = gene_distance 
            min_gene1, min_gene2 = gene1_name, gene2_name
    return min_distance, min_gene1, min_gene2


def check_gene_coverage(sequence_record, length=None, check_for_overlap=True):
    """ Given a GFF sequence record (from BCBio.GFF parser), determine what fraction of its length is covered by genes.
    Length can be passed separately, or determined from the record (Caution: approximate! First-last feature only).
    """
    if length is None:  length = sequence_record.seq._length

    gene_length_total = 0
    for gene in sequence_record.features:
        gene_length_total += gene.location.end.position - gene.location.start.position
    gene_coverage_fraction = float(gene_length_total)/length

    # Check for overlapping genes and print a warning, since overlapping genes will make the measurement inaccurate
    if check_for_overlap:
        if check_for_overlapping_genes(sequence_record):
            print "WARNING: There are overlapping genes! %% of length covered by genes may not be accurate."
    # MAYBE-TODO actually adjust the measurement for overlapping genes?  Nah, too much work, not enough need for now.

    return gene_coverage_fraction 


######### Test code #########

class Testing_(unittest.TestCase):
    """ Unit-tests for _____. """
    def test__(self):
        print "UNIT-TEST NOT IMPLEMENTED"
    # MAYBE-TODO implement - really a run-test may be enough here...


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    print "RUN-TEST NOT IMPLEMENTED"
    # TODO implement!  The input file should be test_data/test_reference.gff3; write expected output!
    # see mutant_count_alignments.py do_test_run for how this can be done


######### Main function code #########

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as usage."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-r', '--record_structure', action="store_true", default=False, 
                      help='Show the record structures (for example gene->mRNA->CDS/UTR). Default %default') 
    parser.add_option('-R', '--no_record_structure', action="store_false", dest='record_structure')

    parser.add_option('-c', '--feature_type_counts', action="store_true", default=True, 
                      help='Count the number of feature types in file (gene, mRNA, exon, etc). Default %default') 
    parser.add_option('-C', '--no_feature_type_counts', action="store_false", dest='feature_type_counts')

    parser.add_option('-g', '--gene_counts', action="store_true", default=False, 
                      help="Count genes per chromosome, and the approximate fraction of each chromosome covered by genes. "
                          +"Default %default") 
    parser.add_option('-G', '--no_gene_counts', action="store_false", dest='gene_counts')
    parser.add_option('-d', '--print_seq_details', action="store_true", default=False, 
                      help='Print full GFF details for each chromosome (only if -g). Default %default') 
    parser.add_option('-D', '--no_print_seq_details', action="store_false", dest='print_seq_details')

    parser.add_option('-o', '--check_gene_overlaps', action="store_true", default=True, 
                      help='Check for overlapping genes, distances, ID uniqueness, etc; count genes. Default %default') 
    parser.add_option('-O', '--no_check_gene_overlaps', action="store_false", dest='check_gene_overlaps')

    parser.add_option('-f', '--gene_feature_structure_counts', action="store_true", default=True, 
                      help='Display gene counts by UTR/exon count, order, etc; check feature overlaps. Default %default') 
    parser.add_option('-F','--no_gene_feature_structure_counts', action="store_false",dest='gene_feature_structure_counts')
    parser.add_option('-u', '--full_feature_structures', action="store_true", default=False, 
                      help='With -f option, show full as well as simplified feature structures. Default %default') 
    parser.add_option('-U','--no_full_feature_structures', action="store_false",dest='full_feature_structures')
    parser.add_option('-n', '--genes_to_display', type="int", default=5, metavar='N', 
                      help="When showing gene counts per group (-f), show N example genes (-1: all) (default %default) ")
    parser.add_option('-e', '--exon_number_cutoff', type="int", default=30, metavar='N', 
                      help="When categorizing genes by exon number, lump together all above N (default %default) ")
    parser.add_option('-Y', '--N_detail_run_groups', type="int", default=5, metavar='N', 
                      help="How many passes to split reading the file into with -f option (default %default) "
                          +"- may take a lot of memory (and CPU) if read in a single pass; too many passes waste CPU.")

    parser.add_option('-s', '--source_counts', action="store_true", default=False, 
                      help='Count the features by source (not very useful unless file is mixed-source). Default %default') 
    parser.add_option('-S', '--no_source_counts', action="store_false", dest='source_counts')

    parser.add_option('-l', '--all_gff_limits', action="store_true", default=False, 
                      help='Output all feature counts: by type, source (-cs), chromosome, maybe other? Default %default')
    parser.add_option('-L', '--no_all_gff_limits', action="store_false", dest='all_gff_limits')

    parser.add_option('-E', '--everything', action='store_true', default=False, 
                      help="Run all tests; ignore other options")

    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on test input file, check output against reference. Ignores all other options/arguments.")

    return parser


# Help function for making sure the eval/exec-based -E implementation is right
_sort_dict_string = lambda dict_repr: sorted(str(dict_repr).strip('{}').split(', '))


def run_main_function(infile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument is generated by an optparse parser.
    """
    examiner = GFF.GFFExaminer()

    # turn all True/False options to True (requires eval magic to detect all options)
    if options.everything:
        option_dict = eval(str(options))
        error_text = "The -E option isn't working right, turn everything on by hand!"
        assert _sort_dict_string(option_dict) == _sort_dict_string(options), error_text
        for (key,value) in option_dict.items():
            if value is False:      # must use 'is', not '==': '0==False' is true, but '0 is False' isn't
                exec('options.%s = True'%key)
                assert eval('options.%s'%key) is True, error_text


    # I'm not sure why I need to open the file separately for each operation, but it doesn't work otherwise...

    if options.record_structure:
        with open(infile) as INFILE:
            print "\n *** Record structures ***"
            pprint.pprint(examiner.parent_child_map(INFILE))

    if options.feature_type_counts:
        with open(infile) as INFILE:
            print "\n *** Type counts ***"
            pprint.pprint({'gff_type': examiner.available_limits(INFILE)['gff_type']})

    if options.source_counts:
        print "\n *** Source and source/type counts ***"
        with open(infile) as INFILE:
            pprint.pprint({'gff_source': examiner.available_limits(INFILE)['gff_source']})
        with open(infile) as INFILE:
            pprint.pprint({'gff_source_type': examiner.available_limits(INFILE)['gff_source_type']})

    if options.all_gff_limits:
        with open(infile) as INFILE:
            print "\n *** All GFF file limit field values ***"
            pprint.pprint(examiner.available_limits(INFILE))

    if options.gene_counts or options.print_seq_details or options.check_gene_overlaps:
        if options.gene_counts or options.print_seq_details:
            print "\n *** Gene and other data per chromosome ***"
        if options.gene_counts:
            print "       (Caution: approximate sequence length is calculated to last gene only!)"
        if options.gene_counts or options.print_seq_details:
            print ""
        if options.check_gene_overlaps:
            total_chromosomes = 0
            total_genes = 0
            overlapping_gene_pairs = []
            min_gene_distance_data = None
            gene_IDs = []
        with open(infile) as INFILE:
            for record in GFF.parse(INFILE, limit_info={'gff_type': ['gene']}):
                if options.gene_counts or options.print_seq_details:
                    print " * sequence %s: %s genes"%(record.id, len(record.features))
                if options.gene_counts:
                    gene_coverage_percent = "%.0f%%"%(100*check_gene_coverage(record,None,False))
                    print "Approximate length %s bp, with %s covered by genes"%(record.seq._length, gene_coverage_percent)
                    # MAYBE-TODO get the chromosome length from genome fasta file - I think GFF parser can add those?
                if options.print_seq_details:               
                    print "GFF parser details:"
                    print record
                if options.gene_counts or options.print_seq_details:
                    print ""
                if options.check_gene_overlaps:
                    total_chromosomes += 1
                    total_genes += len(record.features)
                    overlapping_gene_pairs += check_for_overlapping_genes(record)
                    min_gene_distance_data = find_min_gene_distance(record, min_gene_distance_data)
                    gene_IDs += [gene.id for gene in record.features]

    if options.check_gene_overlaps:
        print "\n *** Gene overlaps ***"
        print "Total %s genes on %s chromosomes."%(total_genes, total_chromosomes)
        for geneA,geneB in overlapping_gene_pairs:  print "Overlapping gene pair!  IDs: %s, %s."%(geneA,geneB)
        if not overlapping_gene_pairs:              print "No overlapping genes."
        print "Minimum distance between two genes is %s (genes %s, %s)."%min_gene_distance_data 
        IDs_are_unique = True
        for gene_ID,count in count_list_values(gene_IDs).iteritems():
            if count>1:    
                print "Non-unique gene ID! Gene %s occurs %s times."%(gene_ID,count)
                IDs_are_unique = False
        if IDs_are_unique is True: 
            print "All gene IDs are unique."


    if options.gene_feature_structure_counts:
        print "\n *** Gene counts by feature structure ***"
        genes_by_feature_structure = defaultdict(set)       # set() returns an empty set

        ## Go over subsets of chromosomes at once, to avoid reading the whole file into memory at once
        # First get the list of all chromosomes in the file, WITHOUT reading it all into memory
        with open(infile) as INFILE:
            GFF_limit_data = examiner.available_limits(INFILE)
            chromosomes_and_counts = dict([(c,n) for ((c,),n) in GFF_limit_data['gff_id'].items()])

        # Now lump the chromosomes into N_run_groups sets with the feature counts balanced between sets, 
        #  to avoid using too much memory (by reading the whole file at once), 
        #   or using too much time (by reading the whole file for each chromosome/scaffold)
        chromosome_sets = split_into_N_sets_by_counts(chromosomes_and_counts, options.N_detail_run_groups)
        overlapping_feature_genes = []

        ### go over all mutants on each chromosome, figure out which gene they're in (if any), keep track of totals
        # keep track of all the mutant and reference chromosomes to catch chromosomes that are absent in reference
        for chromosome_set in chromosome_sets:
            genefile_parsing_limits = {'gff_id': list(chromosome_set)}
            with open(infile) as INFILE:
                for chromosome_record in GFF.parse(INFILE, limit_info=genefile_parsing_limits):
                    for gene in chromosome_record.features:
                        if len(gene.sub_features)==0:
                            genes_by_feature_structure[('NO_mRNA',)].add(gene.id)
                        elif len(gene.sub_features)>1:
                            genes_by_feature_structure[('MULTIPLE_mRNAs',)].add(gene.id)
                        else:
                            [mRNA] = gene.sub_features
                            if check_for_overlapping_features(mRNA, gene.id):
                                overlapping_feature_genes.append(gene.id)
                            if not mRNA.type=='mRNA':
                                genes_by_feature_structure[('NON_mRNA_PRIMARY_FEATURE',)].add(gene.id)
                            elif len(mRNA.sub_features)==0:
                                genes_by_feature_structure[('NO_mRNA_SUBFEATURES',)].add(gene.id)
                            else:
                                # The features are NOT SORTED in a normal gff file!!  Need to sort them.
                                features_by_pos = sorted([(f.location.start.position,f.type) for f in mRNA.sub_features])
                                # reverse the structure if the gene's on the minus strand, since it's sorted by position
                                if mRNA.strand in [-1, '-']:    features_by_pos.reverse()
                                feature_structure = tuple([t for (p,t) in features_by_pos])
                                genes_by_feature_structure[feature_structure].add(gene.id)

        for gene in overlapping_feature_genes:      print "Overlapping features in gene!  Gene ID: %s."%gene
        if not overlapping_feature_genes:           print "No genes have overlapping features."

        # make new dict with all rows of 'CDS' changed to a single 'exon/s' to end up with fewer structure variants, 
        #   also count genes with different exon numbers
        exon_number_gene_counts = defaultdict(set)              # passing the set function: equivalent to lambda: set()
        UTR_5prime_number_gene_counts = defaultdict(set)
        UTR_3prime_number_gene_counts = defaultdict(set)
        genes_by_simpler_feature_structure = defaultdict(set)
        for feature_structure,gene_set in genes_by_feature_structure.iteritems():
            simple_feature_structure = []
            exon_count, UTR_5prime_count, UTR_3prime_count = 0,0,0
            for feature in feature_structure:
                if feature=='CDS':
                    exon_count += 1
                    if simple_feature_structure==[] or not simple_feature_structure[-1]=='exon/s':
                        simple_feature_structure.append('exon/s')
                elif feature=='five_prime_UTR':
                    UTR_5prime_count += 1
                    if simple_feature_structure==[] or not simple_feature_structure[-1]=="5'UTR/s":
                        simple_feature_structure.append("5'UTR/s")
                elif feature=='three_prime_UTR':
                    UTR_3prime_count += 1
                    if simple_feature_structure==[] or not simple_feature_structure[-1]=="3'UTR/s":
                        simple_feature_structure.append("3'UTR/s")
                else:
                    simple_feature_structure.append(feature)
            genes_by_simpler_feature_structure[tuple(simple_feature_structure)] |= gene_set
            UTR_5prime_number_gene_counts[UTR_5prime_count] |= gene_set
            UTR_3prime_number_gene_counts[UTR_3prime_count] |= gene_set
            if exon_count >= options.exon_number_cutoff:
                exon_number_gene_counts['%s+'%options.exon_number_cutoff] |= gene_set
            else:
                exon_number_gene_counts[exon_count] |= gene_set

        N = options.genes_to_display if options.genes_to_display>=0 else None

        print( "\n * Gene counts by simplified feature structure (adjacent exons/UTRs combined) "
              +"(total %s structures)"%len(genes_by_simpler_feature_structure))
        genes_by_simpler_feature_structure = [(len(s),n,s) for (n,s) in genes_by_simpler_feature_structure.items()]
        genes_by_simpler_feature_structure.sort(reverse=True)
        for count,feature_structure,gene_set in genes_by_simpler_feature_structure:
            print "%s genes [%s]:  %s"%(count, ', '.join(feature_structure), ', '.join(list(gene_set)[:N]))

        print "\n * Gene counts by exon number (total %s values)"%len(exon_number_gene_counts)
        for exon_number,gene_set in sorted(exon_number_gene_counts.items()):
            print "%s genes with %s exons:  %s"%(len(gene_set), exon_number, ', '.join(list(gene_set)[:N]))

        print "\n * Gene counts by 5'UTR number (total %s values)"%len(UTR_5prime_number_gene_counts)
        for UTR_number,gene_set in sorted(UTR_5prime_number_gene_counts.items()):
            print "%s genes with %s 5'UTRs:  %s"%(len(gene_set), UTR_number, ', '.join(list(gene_set)[:N]))

        print "\n * Gene counts by 3'UTR number (total %s values)"%len(UTR_3prime_number_gene_counts)
        for UTR_number,gene_set in sorted(UTR_3prime_number_gene_counts.items()):
            print "%s genes with %s 3'UTRs:  %s"%(len(gene_set), UTR_number, ', '.join(list(gene_set)[:N]))

        if options.full_feature_structures:
            genes_by_feature_structure = [(len(s),x,s) for (x,s) in genes_by_feature_structure.items()]
            genes_by_feature_structure.sort(reverse=True)
            print "\n * Gene counts by full feature structure (total %s structures)"%len(genes_by_feature_structure)
            for count,feature_structure,gene_set in genes_by_feature_structure:
                print "%s genes [%s]:  %s"%(count, ', '.join(feature_structure), ', '.join(list(gene_set)[:N]))



if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # if ran with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise parse the arguments and run main function
    try:                
        [infile] = args
    except ValueError:  
        parser.print_help()
        sys.exit("\nError: exactly one input gff file required!")

    run_main_function(infile, options)


