#!/usr/bin/env python2.7
""" Convert seq/alignment data into full mutant data, grouped by internal barcode, with positions/genes/annotation.

Outputs a line-per-mutant file containing position info, gene annotation data (optional), total/perfect read count, number of distinct sequences, and some of the sequences/counts (optional), and a summary of read/mutant counts etc. 
Output files are in simple tab-separated plaintext format, with all header/summary lines starting with #.

Grouping reads into mutants is done by clustered IBs (internal barcodes) - clustering the IBs should be done before this.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: mutant_count_alignments.py [options] outfile """
# TODO update docstring!

# basic library
import sys, os, time
import unittest
import pickle
# other packages
import HTSeq
# my modules
from general_utilities import write_header_data
from testing_utilities import run_functional_tests
import mutant_IB_RISCC_classes

######################################################### Main function code ####################################################

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    ### Input files and basic descriptions
    # MAYBE-TODO make it possible to put multiple files for all these?
    # TODO figure out what options should be used when parsing RISCC data vs IB-only data, and how to merge them together!
    parser.add_option('-c', '--casette_side_reads', default=None, metavar='FILE', 
                      help="SAM file containing aligned cassette-side reads (from RISCC or ChlaMmeSeq protocol or such; "
                          +"orientation defined with -d, cassette end defined with -e). Required.")
    parser.add_option('-b', '--internal_barcode_reads', default=None, metavar='FILE', 
                      help="Fastq file containing internal barcode reads (from RISCC or IB-only PCR protocol), "
                          +"(read IDs should be the same as matching cassette-side reads). Required.")
    parser.add_option('-B', '--internal_barcode_clusters', default=None, metavar='FILE', 
                      help="Python pickle file containing internal barcode clustering data (from RISCC or IB-only PCR protocol), "
                          +"as a centroid_seq:set_of_IB_seqs dictionary (generated by ____). Required.")
    parser.add_option('-g', '--genome_side_reads', default=None, metavar='FILE', 
                      help="SAM file containing aligned genome-side reads (from RISCC protocol or such; "
                          +"in orientation toward the cassette). "
                          +"Should have same read IDs as matching cassette-side reads, and be from the same side. Optional.")
    # MAYBE-TODO make the IB file optional and allow the cluster file to be based on cassette-side seqs? Is that a good idea?
    # MAYBE-TODO make the IB+cluster files optional, and use the old clustering-by-alignment method in that case?
    parser.add_option('-e', '--read_cassette_end', choices=mutant_IB_RISCC_classes.SEQ_ENDS, default='5prime', 
                      metavar='|'.join(mutant_IB_RISCC_classes.SEQ_ENDS), 
                      help="Which end of the cassette are the sequenced reads from? (default %default).")
    parser.add_option('-d','--relative_read_direction', choices=mutant_IB_RISCC_classes.RELATIVE_READ_DIRECTIONS, default='outward',
                      metavar='|'.join(mutant_IB_RISCC_classes.RELATIVE_READ_DIRECTIONS), 
                      help="Are the cassette-side reads oriented inward or outward to the cassette? (default %default).")

    ### Functionality options
    parser.add_option('-X', '--best_genome_side_only', action="store_true", default=False,
                      help="Instead of storing all genome-side reads per mutant, only store the 'best' one (uniquely aligned, "
                          +"to the same general location as cassette-side but with the maximal distance from it) "
                          +"and the overall counts of genome-side read categories, to save memory. (Default %default).")
    # gene-finding and gene-annotation options 
    parser.add_option('-A', '--gene_annotation_folder', default='None', metavar='FOLDER', 
                      help="Folder containing gene position and annotation files (files should be from Phytozome: "
                          +"gff file with gene positions, *_annotation_info.txt, maybe *_defline.txt etc) "
                          +"if None, gene IDs and annotation won't be added (default %default)")
    # MAYBE-TODO change default to be based on an environmental variable?
    # TODO or would it be better to still give the gene/annotation file names, but ALSO check them against the genome version?
    # MAYBE-TODO add extra annotation files of some sort? Martin had lots of ideas...
    parser.add_option('-G', '--genome_version', type='int', default=5, metavar='G', 
                      help="Which genome version the input files were aligned against - picks the matching gene position "
                          +"and annotation/etc files from the -A folder (4 for v4.* genome, 5 for v5.*, etc) (default %default)")
    # TODO pick the correct gene/annotation files out of the folder; make sure it's consistent with the alignment genome!
    parser.add_option('--detailed_gene_features', action="store_true", default=True,
                      help="Find out what part of the gene (UTR,intron,exon) a mutant hit, based on the -g file "
                          +"(default %default). May take a lot of memory - increase --N_detail_run_groups to fix that.")
    parser.add_option('--no_detailed_gene_features', action="store_false", dest='detailed_gene_features',
                      help="Turns --detailed_gene_features off.")
    parser.add_option('--N_detail_run_groups', type="int", default=5, metavar='N', 
                      help="How many passes to split reading the detailed_gene_features into (default %default) "
                          +"- may take a lot of memory (and CPU) if read in a single pass; too many passes waste CPU.")
    # MAYBE-TODO do we ever want to use --no_detailed_gene_features, really?
    # MAYBE-TODO add a "flank" option (with variable size), to catch mutants that are in the flanks of genes? Do we care?
    # MAYBE-TODO add "distance from gene start" and "distance from gene end" fields.

    ### output format options
    parser.add_option('-o', '--sort_data_key', choices=['position','read_count','none'], default='position', 
                      metavar='position|read_count|none', help="Sort the output data: by alignment position, read count, "
                         +"or don't sort at all (default %default) - sorting may be slow for large datasets!")

    parser.add_option('-V', '--verbosity_level', action="store_true", default=1, 
                      help="How much information to print to STDOUT: 0 - nothing, 1 - summary only, "
                          +"2 - summary and progress reports. (Default %default).")
    parser.add_option('-q', '--quiet', action="store_const", const=0, dest='verbosity_level', help="Equivalent to -V 0.")
    parser.add_option('-v', '--verbose', action="store_const", const=2, dest='verbosity_level', help="Equivalent to -V 2.")

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    return parser

# TODO merge this with a similar code bit in mutant_join_datasets.py to minimize code duplication
def save_dataset_files(dataset, outfile, verbosity_level=0, if_pickle=True, count_cassette=True, count_other=True, 
                       sort_data_by='position', options="N/A"):
    """ Print summary and data to output file; optionally print summary to stdout; optionally pickle dataset to picklefile. 
    
    The options argument is only used to be printed in the header to make it clear how the file was generated - 
     it should be the applicable optparse options object if there is one, or a text message otherwise.
    """
    # print summary info to stdout if desired
    if verbosity_level>1: print "\nDATA SUMMARY:"
    if verbosity_level>0: dataset.print_summary(count_cassette=count_cassette, count_other=count_other)
    # print full data to outfile
    if verbosity_level>1: print "printing output - time %s."%time.ctime()
    with open(outfile,'w') as OUTFILE:
        write_header_data(OUTFILE,options)
        OUTFILE.write("### SUMMARY:\n")
        dataset.print_summary(OUTFILE, line_prefix="#  ", header_prefix="## ", 
                              count_cassette = count_cassette, count_other=count_other)
        OUTFILE.write("### HEADER AND DATA:\n")
        dataset.print_data(OUTPUT=OUTFILE, sort_data_by=sort_data_by, header_line=True, header_prefix='# ')
    # print pickled dataset to picklefile, if desired
    if if_pickle:
        outfile_basename = os.path.splitext(outfile)[0]
        pickled_outfile = outfile_basename + '.pickle'
        with open(pickled_outfile,'w') as PICKLEFILE:
            pickle.dump(dataset, PICKLEFILE, 0)


def main(outfile, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """
    outfile_basename = os.path.splitext(outfile)[0]
    # MAYBE-TODO change outfile to a folder?  Since there are 3 output files at this point...

    ### generate empty alignment set object with basic read position/orientation properties defined by options
    all_alignment_data = mutant_IB_RISCC_classes.Insertional_mutant_pool_dataset(options.read_cassette_end, 
                                                                                 options.relative_read_direction)
    ### parse input files and store data - the add_RISCC_alignment_files_to_data function here does pretty much all the work!
    all_alignment_data.add_RISCC_alignment_files_to_data(options.internal_barcode_clusters, options.internal_barcode_reads, 
                                                         options.casette_side_reads, options.genome_side_reads, 
                                                         best_genome_side_only = options.best_genome_side_only)

    ### MAYBE-TODO some kind of metadata parsing?  If so, copy relevant options and code from mutant_count_alignments.py.
    # TODO do we want to deal with metadata in the new IB+Carette pipeline? How?
    # This is mostly relevant to the cassette-side file, since that one had some reads removed due to not matching the expected 
    #  structure, etc.  But there can also be metadata for the IB clustering, and for genome-side alignment... 
    #  MAYBE-TODO include that in the same file or a separate one? Do we even care?

    # I'm assuming input isn't collapsed to unique - it makes no sense with RISCC or IB data, since read IDs are important! 
    # MAYBE-TODO unless I write something where we can do the alignments on collapsed-to-unique seqs but the look up the read IDs
    #  by sequence?  Probably not very useful, since the genome-side seqs will mostly all be unique anyway, 
    #  and doing this only to the cassette-side reads would save less than half of the alignment work.

    ### optionally merge some mutant categories?  Or at least count them?
    # It makes no sense to merge adjacent mutants when the mutants are based on IB clustering rather than position.
    # MAYBE-TODO implement merging opposite-tandem mutants (two sides of same two-cassette insertion) based on position?  
    #  That means merging two different IBs... This should be done for same-side cases and ALSO for 3'+5' cases, 
    #  which we don't have any of yet, so can probably skip it for now...
    # If we do merging, change the counting below to print the data to MERGEFILE/CASSETTE_MERGEFILE rather than sys.stdout!

    ### optionally parse gene position/info/annotation files and look up the genes for each mutant in the data
    if options.gene_annotation_folder not in (None, 'None', 'NONE', 'none', ''):
        # TODO find gene gff3 file and any annotation files based on options.gene_annotation_folder and options.genome_version
        if options.verbosity_level>1: print "adding genes from file %s to mutant data - time %s."%(genefile, time.ctime())
        all_alignment_data.find_genes_for_mutants(genefile, detailed_features=options.detailed_gene_features, 
                                                  N_run_groups=options.N_detail_run_groups, verbosity_level=options.verbosity_level)
        if options.verbosity_level>1: 
            print "adding gene annotation from file %s - time %s."%(options.gene_annotation_file, time.ctime())
        all_alignment_data.add_gene_annotation(gene_annotation_file, 
                   if_standard_Phytozome_file=options.annotation_file_standard_type, print_info=(options.verbosity_level >= 2))

    ### output data to files
    save_dataset_files(all_alignment_data, outfile, options.verbosity_level, True, True, True, 
                       True, options.sort_data_key, options)
    # TODO write some info about all the other files that go with this one (pickle, merging-info, *cassette*)


# TODO need between-mutant sanity checks!  Like cases with the same cassette-side seq in different mutants (when it's genome-aligned), and cases with very close positions, and such.

# TODO write separate program for removing some mutants based on other file


######################################################### Testing and main ####################################################

# TODO make new tests for new IB+Carette pipeline!
def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    test_folder = "test_data"
    aln_infile0 = "test_data/INPUT_alignment0_old-format.sam"
    aln_infile1 = "test_data/INPUT_alignment1_genomic-unique.sam"
    aln_infile2 = "test_data/INPUT_alignment2_for-genes.sam"
    gff_genefile = "test_data/INPUT_gene-data-1_all-cases.gff3"
    dataset_to_remove = "test_data/INPUT_mutants_to_remove.txt"

    test_runs = [
                 ('cassette-end-5prime', "-e 5prime -r forward -n3 -L", [aln_infile1]),
                 ('cassette-end-3prime', "-e 3prime -r forward -n3 -L", [aln_infile1]),
                 ('read-direction-reverse', "-r reverse -e 5prime -n3 -L", [aln_infile1]),
                 ('sorted-by-count', "-o read_count -e 5prime -r forward -n3 -L", [aln_infile1]),
                 ('with-gene-info_merged', "-e 5prime -r forward -g %s -n0"%gff_genefile, [aln_infile2]),
                 ('remove-from-other-all', "-x %s -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-from-other-min4', "-x %s -z4 -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-from-other-perfect', "-x %s -p -z4 -n0"%dataset_to_remove, [aln_infile2]),
                 ('remove-not-other-all', "-X %s -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-not-other-min4', "-X %s -Z4 -n0"%dataset_to_remove, [aln_infile2]), 
                 ('remove-not-other-perfect', "-X %s -P -Z4 -n0"%dataset_to_remove, [aln_infile2]),
                ]
    # TODO remove the mutation-detection lines from aln_infile1 and from all the outfiles, since that's been simplified?
    # TODO add run-test for removing data from multiple files?
    # MAYBE-TODO add run-test for a metadata file with 5' and 3' read counts?
    # MAYBE-TODO add run-tests for other mutant-merging options?  But they all have pretty good unit-tests.
    # MAYBE-TODO add run-test for --gene_annotation_file?
    # MAYBE-TODO add run-test for --input_collapsed_to_unique?  Or is that unit-tested already?

    # convert tests into (testname, arg_and_infile_string) format, adding the options that are always used
    test_names_and_args = [('count-aln__'+testname, test_args+' -q '+' '.join(infiles)) 
                           for testname,test_args,infiles in test_runs]

    parser = define_option_parser()
    argument_converter = lambda parser,options,args: (args[:-1], args[-1], options)
    return run_functional_tests(test_names_and_args, parser, main, test_folder, 
                                argument_converter=argument_converter, append_to_outfilenames='.txt') 


class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__(self):
        print "NO UNIT-TESTS FOR THIS MODULE"


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    parser = define_option_parser()
    (options, args) = parser.parse_args()

    # if run with -t or -T option, do the relevant tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        print("\n * unit-tests for the mutant_IB_RISCC_classes.py module")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(mutant_IB_RISCC_classes)
        unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        print("\n * unit-tests for this module (%s)"%sys.argv[0])
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()
    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise parse the arguments and run main function
    try:
        [outfile] = [args]
    except ValueError:
        parser.print_help()
        sys.exit("\nError: exactly one outfile is required! Infiles should be given as options.")
    main(outfile, options)
