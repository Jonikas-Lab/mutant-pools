#!/usr/bin/env python

import os

def parse_gene_annotation_file(gene_annotation_filename, standard_Cre_file=False, header_fields=None, 
                               genes_start_with='Cre', strip_last_gene_fields=0, gene_field_delimiter='.', 
                               pad_with_empty_fields=True, ignore_comments=False, verbosity_level=1):
    """ Parse tab-separated gene annotation file; return a geneID:data_list dict, and a column header list (if present).
    Try to automatically detect if the file has a header; optionally use header provided; 
     return None for header if none was provided or detected. 
    Remove columns with no contents.
    """
    # TODO update/rewrite docstring!

    if not os.path.lexists(gene_annotation_filename):
        raise Exception("Couldn't find the %s gene annotation file!"%gene_annotation_filename)
    if verbosity_level>0:
        print " *** Parsing file %s for gene annotation info..."%gene_annotation_filename

    # special case for the Creinhardtii_169_annotation_info.txt file and similar (from README for that file)
    if standard_Cre_file:
        header_fields = ['Phytozome transcript name', 'PFAM', 'Panther', 'KOG', 'KEGG ec', 'KEGG Orthology', 
                         'best arabidopsis TAIR10 hit name', 'best arabidopsis TAIR10 hit symbol', 
                         'best arabidopsis TAIR10 hit defline', 
                         'best rice hit name', 'best rice hit symbol', 'best rice hit defline']
        strip_last_gene_fields = 2

    # parse the whole file into lists of tab-separated fields
    # MAYBE-TODO this reads in the whole file at once - not very efficient! Change to a generator?  
    #  (except we need special treatment for the first line which may or may not be a header line...)
    data_by_row = []
    for line in open(gene_annotation_filename):
        if ignore_comments and line[0]=='#':    continue
        fields = line.strip().split('\t')
        data_by_row.append(fields)
        
    # if the first line doesn't start with an expected gene ID, assume it's a header
    if not data_by_row[0][0].startswith(genes_start_with):
        header = data_by_row[0]
        del data_by_row[0]
        if verbosity_level>0:
            print "Assuming the first line is a header: %s"%'\t'.join(header)
    else: 
        header = None
    # if any of the other lines doesn't start with a Cre* gene ID, fail!
    for row in data_by_row:
        if not row[0].startswith(genes_start_with):
            raise Exception("Can't parse file %s - found non-header line that doesn't start "%gene_annotation_filename
                            +"with a %s gene ID!\n  \"%s\""%(genes_start_with, '\t'.join(row)))
        # TODO do I want to somehow prevent silently replacing file header with header from arguments?...
    # use the header from arguments if provided
    if header_fields is not None:
        header = header_fields
        if verbosity_level>1:
            print "Using header provided: %s"%'\t'.join(header)

    # check that all the data lengths line up (if they don't, don't throw an error
    mismatched_lengths = False
    data_lengths = set([len(row) for row in data_by_row])
    if not len(data_lengths)==1:
        if verbosity_level>1 or (not pad_with_empty_fields and verbosity_level>0):
            print "Warning: not all data rows have the same length! Lengths found: %s"%list(data_lengths)
        mismatched_lengths = True
    if header is not None and len(header) not in data_lengths:
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
    if not mismatched_lengths:
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
        gene = row[0]
        if strip_last_gene_fields:
            gene = gene_field_delimiter.join(gene.split(gene_field_delimiter)[:-strip_last_gene_fields])
        data = row[1:]
        if gene in data_by_gene:
            if verbosity_level>0:
                print "Warning: gene %s appears twice in the data! Using last appearance."%gene
        data_by_gene[gene] = data

    # remove the first word from the header, since it should be "gene ID" or such; 
    #  change spaces to underscores in header fields for readability
    if header:  
        del header[0]
        header = [s.replace(' ','_') for s in header]

    if verbosity_level>0:
        print " *** DONE Parsing gene annotation file"
    return data_by_gene, header

# TODO add unit-tests or something?  Right now I don't really know if this works at all!

