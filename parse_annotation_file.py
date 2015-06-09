#!/usr/bin/env python2.7

# standard library
from __future__ import division
import sys
import unittest
import os
from collections import defaultdict

# Gene annotation file data for v5.5 genome: 
#   list of (filename, content_headers, ID_column, content_columns, if_join_all_later_fields) tuples:
DEFAULT_GENE_ANNOTATION_FILES_v5p5 = [
    ('Creinhardtii_281_v5.5.geneName.txt', ['gene_name'], 0, [1], False),
    ('Creinhardtii_281_v5.5.defline.txt', ['defline'], 0, [1], False),
    ('Creinhardtii_281_v5.5.description.txt', ['description'], 0, [1], False),
    ('Creinhardtii_281_v5.5.synonym.txt', ['synonyms'], 0, [1], True),
    ('Creinhardtii_281_v5.5.annotation_info.txt', 'PFAM Panther KOG KEGG_ec KEGG_Orthology Gene_Ontology_terms best_arabidopsis_TAIR10_hit_name best_arabidopsis_TAIR10_hit_symbol best_arabidopsis_TAIR10_hit_defline'.split(), 
     1, [4, 5, 6, 7, 8, 9, 10, 11, 12], False),
]
DEFAULT_GENE_ANNOTATION_FOLDER = os.path.expanduser('~/experiments/reference_data/chlamy_annotation')
DEFAULT_GENE_ANNOTATION_FILES_v5p5 = [(os.path.join(DEFAULT_GENE_ANNOTATION_FOLDER, f), h, i, c, j) 
                                      for (f, h, i, c, j) in DEFAULT_GENE_ANNOTATION_FILES_v5p5]

# LATER-TODO for synonyms files, strip the "t.*" part from the synonym ID columns as well as the gene ID column?

# Older deprecated annotation info
HEADER_FIELDS_v4 = ['Phytozome transcript name', 
                    'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 
                    'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 'best arabidopsis TAIR10 hit defline', 
                    'best rice hit name', 'best rice hit symbol', 'best rice hit defline']

HEADER_FIELDS_v5 = ['Phytozome internal transcript ID', 
                    'Phytozome gene locus name', 'Phytozome transcript name', 'Phytozome protein name', 
                    'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 'Gene Ontology terms', 
                    'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 'best arabidopsis TAIR10 hit defline', 
                    'best rice hit name', 'best rice hit symbol', 'best rice hit defline']


def parse_gene_annotation_file(gene_annotation_filename, content_header_fields, gene_ID_column=0, content_columns=[1], 
                               if_join_all_later_fields=False, pad_with_empty_fields=True,
                               strip_gene_fields_start=".t", genes_start_with=None, ignore_comments=False, verbosity_level=1):
    """ Parse tab-separated gene annotation file; return gene:annotation_list dictionary.

    Use column gene_ID_column to determine gene IDs; optionally shorten the gene name by truncating it starting with the 
     strip_gene_fields_start value if found (if not None) - e.g. if value is '.t', Cre01.g123450.t2.1 would become Cre01.g123450.
    If genes_start_with is not None, make sure all gene IDs start with it.

    If pad_with_empty_fields is True, pad shorter lines to the max length. 
    If ignore_comments is True, skip lines starting with #.

    Use content_columns for the values, or if if_join_all_later_fields is True, then use a comma-separated list of ALL
     fields starting with the content_columns one.

    Print some info/warnings to stdout depending on verbosity_level (0 - nothing, 1 - some, 2 - max).
    """
    if not os.path.lexists(gene_annotation_filename):
        raise Exception("Couldn't find the %s gene annotation file!"%gene_annotation_filename)
    if verbosity_level>0:
        print "  Parsing file %s for gene annotation info..."%os.path.basename(gene_annotation_filename)

    ### Parse the whole file into lists of tab-separated fields
    #    (could change to a generator, but most of the data has to stay in memory anyway in a different format, so probably no point)
    #    (and we need special treatment for the first line which may or may not be a header line...)
    data_by_row = []
    for line in open(gene_annotation_filename):
        if ignore_comments and line[0]=='#':    continue
        fields = line.strip().split('\t')
        data_by_row.append(fields)
    if verbosity_level>0:
        print "  Parsed %s lines"%len(data_by_row)
        
    # if any of the other lines doesn't start with a Cre* gene ID, fail!
    if genes_start_with is not None:
        for row in data_by_row:
            if not row[gene_ID_column].startswith(genes_start_with):
                raise Exception("Can't parse file %s - found line that doesn't start "%gene_annotation_filename
                                +"with a %s gene ID!\n  \"%s\""%(genes_start_with, '\t'.join(row)))

    # check that all the data lengths line up (if they don't, don't throw an error
    data_lengths = set([len(row) for row in data_by_row])
    if len(data_lengths)>1 and not if_join_all_later_fields:
        mismatched_lengths = True
        if verbosity_level:     print "Not all data rows have the same length! Lengths found: %s"%list(data_lengths)
        if pad_with_empty_fields:
            max_length = max(data_lengths)
            if verbosity_level>0:
                print "Data field numbers vary between rows - padding all lower-length data rows to length %s"%max_length
            for row in data_by_row:
                if len(row)<max_length:
                    row += ['' for x in range(max_length-len(row))]
            mismatched_lengths = False
    else:   mismatched_lengths = False

    # LATER-TODO figure out the total #fields over the whole file, and then: 
    #   - if content_header_fields is None, just use the filename and N empties
    #   - if content_columns is None, use all except the ID column, I guess?
    if len(content_header_fields) != len(content_columns):
        raise Exception("Error: content_header_fields has a different number of fields than content_columns!")
    if max(data_lengths) < max(content_columns):
        raise Exception("content_columns specifies a column that doesn't exist!")

    # MAYBE-TODO remove empty columns (only if all the data lengths match!)

    ### Convert the list-format data into a by-gene dictionary, grabbing only the fields we want
    # we frequently get multiple lines, for different transcripts!  Just concatenate all of them.
    data_by_gene = defaultdict(lambda: [set() for _ in content_columns])
    for data in data_by_row:
        gene = data[gene_ID_column]
        if strip_gene_fields_start is not None:
            gene = gene.split(strip_gene_fields_start)[0]
        if if_join_all_later_fields:
            if len(content_columns) > 1:
                raise Exception("if_join_all_later_fields not implemented with multple content_columns!")
            data_by_gene[gene][0].add(','.join(str(x) for x in data[content_columns[0]:] if x.strip()))
        else:
            for (new_col, old_col) in enumerate(content_columns):
                data_by_gene[gene][new_col].add(data[old_col])

    # At the end, change the sets to strings, and also remove empty strings
    for gene,data in data_by_gene.items():
        data = [', '.join([f for f in fields if f.strip()]) for fields in data]
        data_by_gene[gene] = [x if x else '-' for x in data]

    if verbosity_level>0:
        print "  DONE Parsing gene annotation file - found %s genes"%len(data_by_gene)
    return defaultdict(lambda: ['-' for _ in content_columns], data_by_gene)


def get_all_gene_annotation(genome_version=None, gene_annotation_files=None, print_info=False):
    """ Grab all the annotation (depends on genome version); return gene:annotation_list dict and header list.

    Can provide a gene_annotation_files dict instead - in the same format as DEFAULT_GENE_ANNOTATION_FILES_v5p5 here.
    """
    if genome_version is None and gene_annotation_files is None:
        raise Exception("User has to provide genome_version or gene_annotation_files!")
    if genome_version is not None and gene_annotation_files is not None:
        raise Exception("User should not provide both genome_version and gene_annotation_files, just one!")
    elif gene_annotation_files is None:
        if genome_version == 5.5:
            gene_annotation_files = DEFAULT_GENE_ANNOTATION_FILES_v5p5
        else:
            raise Exception("Genome version %s not implemented right now!"%genome_version)
    # MAYBE-TODO add options for strip_gene_fields_start etc
    gene_annotation_dicts = [parse_gene_annotation_file(filename, content_header_fields, ID_column, content_columns, 
                                                        if_join_all_later_fields, pad_with_empty_fields=True,
                                                        strip_gene_fields_start=".t", genes_start_with=None, 
                                                        ignore_comments=False, verbosity_level=print_info) 
                             for (filename, content_header_fields, ID_column, content_columns, if_join_all_later_fields) 
                             in gene_annotation_files]
    full_header = sum([content_headers for (_, content_headers, _, _, _) in gene_annotation_files], [])

    all_gene_IDs = set.union(*[set(d.keys()) for d in gene_annotation_dicts])
    full_annotation_dict = defaultdict(lambda: ['-' for _ in full_header])
    for gene in all_gene_IDs:
        full_annotation_dict[gene] = sum([d[gene] for d in gene_annotation_dicts], [])
    return full_annotation_dict, full_header


### Bit of old code to get gene names from gff3 file (no longer necessary with v5.5 genome, which has a tab-sep geneName file):
#   new_fields.append('transcript_names')
#   genename_dict = defaultdict(set)
#   for line in open(gff_file_for_gene_names):
#       if line.startswith('#'):    continue
#       all_fields = line.strip().split('\t')
#       if all_fields[2] != 'mRNA': continue
#       data = all_fields[8]
#       fields = dict(x.split('=') for x in data.split(';'))
#       gene = fields['Name'].split('.t')[0]
#       try:                genename_dict[gene].add(fields['geneName'])
#       except KeyError:    pass
#   genename_dict = defaultdict(set, {key:','.join(vals) for (key,vals) in genename_dict.items()})


class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__all(self):
        # do it twice just to make sure the results are the same - earlier I had a bug where they weren't!
        for i in (1,2,3,4):
            if i < 3:
                data, header = get_all_gene_annotation(5.5, print_info=False)
            else:
                data, header = get_all_gene_annotation(gene_annotation_files=DEFAULT_GENE_ANNOTATION_FILES_v5p5, print_info=False)
            self.assertEquals(header, 'gene_name defline description synonyms PFAM Panther KOG KEGG_ec KEGG_Orthology Gene_Ontology_terms best_arabidopsis_TAIR10_hit_name best_arabidopsis_TAIR10_hit_symbol best_arabidopsis_TAIR10_hit_defline'.split())
            self.assertEquals(len(data), 17741)
            # a gene that doesn't exist
            self.assertEquals(data['Cre01.g000000'], ['-' for x in header])
            # a gene that exists
            self.assertEquals(data['Cre01.g000150'], ['ZRT2', 'Zinc-nutrition responsive permease transporter', 'Zinc/iron permease; related to plant homologs; ZIP family, subfamily I; appears as CrZIP2 in PMID: 15710683; PMID: 16766055', 'g6.t1,Cre01.g000150.t1.1', 'PF02535', 'PTHR11040,PTHR11040:SF30', 'KOG1558', '-', '-', 'GO:0016020,GO:0030001,GO:0046873,GO:0055085', 'AT2G04032.1', 'ZIP7', 'zinc transporter 7 precursor'])
            # a gene that exists and has multiple annotation/synonym lines
            self.assertEquals(data['Cre17.g733850'], 
                              ['-', '-', '-', 'g17740.t1,Cre17.g733850.t1.1, g17740.t2', 'PF01391'] + ['-' for x in range(8)])
    # MAYBE-TODO add a v4 unit-test too, and for other formats?  Once I have them implemented.


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
