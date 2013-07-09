#!/usr/bin/env python

import os

def parse_gene_annotation_file(gene_annotation_filename, standard_Phytozome_file=False, header_fields=None, gene_ID_column=0, 
                               genes_start_with=None, remove_empty_columns=True, strip_gene_fields_start=None,
                               pad_with_empty_fields=True, ignore_comments=False, verbosity_level=1):
    """ Parse tab-separated gene annotation file; return a geneID:data_list dict, and a column header list (if present).

    Try to automatically detect if the file has a header; optionally use header provided; 
     return None for header if none was provided or detected. 
    
    If standard_Phytozome_file is non-zero, use special case headers for standard Phytozome files:
        v4 annotation file (Creinhardtii_169_annotation_info.txt) if value is 4, v5 (Creinhardtii_236_readme.txt) if 5.
        (those files don't have included headers! And the formats sometimes change, so this isn't a guarantee.)

    Optionally shorten the gene name by truncating it starting with the strip_gene_fields_start value if found (if not None):
        for example if strip_gene_fields_start is '.t', Cre01.g123450.t2.1 would become Cre01.g123450.

    Optionally remove empty columns (i.e. fields that are empty in all lines).
    """
    if not os.path.lexists(gene_annotation_filename):
        raise Exception("Couldn't find the %s gene annotation file!"%gene_annotation_filename)
    if verbosity_level>0:
        print " *** Parsing file %s for gene annotation info..."%gene_annotation_filename

    # special cases for the standard Phytozome annotation files
    if standard_Phytozome_file:
        strip_gene_fields_start = '.t'
        if standard_Phytozome_file==4:
            gene_ID_column = 0
            header_fields = ['Phytozome transcript name', 'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 
                             'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 
                             'best arabidopsis TAIR10 hit defline', 
                             'best rice hit name', 'best rice hit symbol', 'best rice hit defline']
        elif standard_Phytozome_file==5:
            gene_ID_column = 1
            header_fields = ['Phytozome transcript name', 
                             'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 'Gene Ontology terms', 
                             'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 
                             'best arabidopsis TAIR10 hit defline', 
                             'best rice hit name', 'best rice hit symbol', 'best rice hit defline']
        else:
            raise Exception("Invalid value for standard_Phytozome_file type arg! %s given, 4/5 accepted."%standard_Phytozome_file)

    # parse the whole file into lists of tab-separated fields
    # MAYBE-TODO this reads in the whole file at once - not very efficient! Change to a generator?  
    #  (except we need special treatment for the first line which may or may not be a header line...)
    data_by_row = []
    for line in open(gene_annotation_filename):
        if ignore_comments and line[0]=='#':    continue
        fields = line.strip().split('\t')
        data_by_row.append(fields)
        
    # if header not given, assume the first line is a header
    if header_fields is None:
        header_fields = data_by_row[0]
        del data_by_row[0]
        if verbosity_level>0:
            print "Assuming the first line is a header: %s"%'\t'.join(header_fields)
    # if any of the other lines doesn't start with a Cre* gene ID, fail!
    if genes_start_with is not None:
        for row in data_by_row:
            if not row[0].startswith(genes_start_with):
                raise Exception("Can't parse file %s - found line that doesn't start "%gene_annotation_filename
                                +"with a %s gene ID!\n  \"%s\""%(genes_start_with, '\t'.join(row)))

    # check that all the data lengths line up (if they don't, don't throw an error
    mismatched_lengths = False
    data_lengths = set([len(row) for row in data_by_row])
    if not len(data_lengths)==1:
        if verbosity_level>1 or (not pad_with_empty_fields and verbosity_level>0):
            print "Warning: not all data rows have the same length! Lengths found: %s"%list(data_lengths)
        mismatched_lengths = True
    if len(header_fields) not in data_lengths:
        if verbosity_level>1 or (not pad_with_empty_fields and verbosity_level>0):
            print("Warning: header has a different number of fields than the data! "
                  +"Header length: %s. Data lengths: %s"%(len(header),list(data_lengths)))
        mismatched_lengths = True
    if len(data_lengths)>1 and pad_with_empty_fields:
        max_length = max(max(data_lengths), len(header))
        if verbosity_level>0:
            print "Data field numbers vary between rows - padding all lower-length data rows to length %s"%max_length
        for row in data_by_row:
            if len(row)<max_length:
                row += ['' for x in range(max_length-len(row))]
        if len(header)<max_length:
            header += ['?' for x in range(max_length-len(header))]
        mismatched_lengths = False

    # remove empty columns (only if all the data lengths match!)
    if remove_empty_columns and not mismatched_lengths:
        data_length = len(header)
        columns_to_remove = []
        for pos in range(data_length):
            values = set([row[pos] for row in data_by_row])
            if len(values)==1:
                value = values.pop()
                if value.strip()=='':
                    if verbosity_level>0:
                        if header:  print "Column %s (%s) is always empty - removing it."%(pos+1, header[pos])
                        else:       print "Column %s is always empty - removing it."%(pos+1)
                    columns_to_remove.append(pos)
                else:
                    if verbosity_level>0:
                        print "Note: all the values in column %s are the same! (%s)"%(pos+1, value)
        for pos in sorted(columns_to_remove, reverse=True):
            if header:  del header[pos]
            for row in data_by_row: 
                del row[pos]

    # convert the list-format data into a by-gene dictionary
    data_by_gene = {}
    for row in data_by_row:
        gene = row[gene_ID_column]
        data = row[0:gene_ID_column] + row[gene_ID_column+1:]
        if strip_gene_fields_start is not None:
            gene = gene.split(strip_gene_fields_start)[0]
        if gene in data_by_gene:
            if verbosity_level>0:
                print "Warning: gene %s appears twice in the data! Using the later appearance."%gene
        data_by_gene[gene] = data

    # remove the first word from the header, since it should be "gene ID" or such; 
    #  change spaces to underscores in header fields for readability
    if header:  
        del header[gene_ID_column]
        header = [s.replace(' ','_') for s in header]

    if verbosity_level>0:
        print " *** DONE Parsing gene annotation file"
    return data_by_gene, header

# TODO add unit-tests or run/tests?  It's not exactly complicated...
