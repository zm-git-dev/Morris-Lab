import csv


# Read the platform map
# 1	DV038825	243377
# 3	BM877548	76748
# 6	BC057190	104625
#
# The platform file allows us to map from microarray spot
# identifier (a small integer) to a variety of gene annotations,
# including various gene identifiers (unigene, refgene, etc) and
# text annotation ("spermin binding protein")
#
# By default we map each id to Entrez GENE ID.
# this behavior can be overridden by passing an optional argument to this routine.
#

def read_platform(reference, column="GENE"): 
    platform = dict()
    genemap = csv.DictReader(open(reference), delimiter='\t',
                             fieldnames=[
                                 'ID', 'MetaCol', 'MetaRow', 'Column', 'Row', 'Biosequence Type',
                                 'Strand Type', 'Region', 'GB_ACC', 'SPOT_ID', 'UNIGENE',
                                 'Reporter Usage', 'Control Type', 'MPEDB Spot ID', 'Designation',
                                 'Related Gene Symbol', 'Description', 'GENE', 'Chromosome',
                                 'Name'])
    for row in genemap:
        try:
            platform[int(row['ID'])] = row['GENE']
        except ValueError:
            continue

    return platform



#
# Read sequencing data
# 
# This should be the output of cufflinks run against the refseq knowgene database.
#
# gene_id	bundle_id	chr	left	right	FPKM	FPKM_conf_lo	FPKM_conf_hi	status
def read_seq_list(seq_file):

    seq_list = list()
    seq_reader = csv.DictReader(open(seq_file), delimiter='\t')

    for row in seq_reader:
        try:
            if (float(row['FPKM']) > 0):
                seq_list.append([row['gene_id'], float(row['FPKM'])])
            continue
        except ValueError:
            continue

    return seq_list


