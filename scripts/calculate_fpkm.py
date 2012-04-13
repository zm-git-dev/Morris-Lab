#! /usr/bin/python
# calculate_fpkm.py
#
# Calculate FPKM values for each gene to which a set of reads has been mapped.
# This script reads from its standard input and sends results to its standard output.
#
# usage: python ${SCRIPTS}/calculate_fpkm.py < tophat_out/aligned_position_stats.txt
# options:
#	-n <depth>	minimum average depth.  Ignore any genes that do not have at
#			least <depth> mapped reads.   Note that mindepth
#			is different from minimum depth;  Depth is measured at
#			a single locus.   depth is measured across the entire
#			gene.
#
# The input file for this script typically comes from rpos_dist.py
#
# input format looks like this:
#    read name       refseq              gene     genelen    region    %    gene%
# ILLUMINA-D43DB5  NM_001164655    9530053A07Rik  3083        CDS    39.80  39.16
#
# This script is only interested in three of these columns:
# 	read_name -- unique name assigned to thie read by the illumina sequencer
# 	refseq	-- refseq name for the gene that this read maps to.
# 	genelen -- length of this gene, in bases (including 5'UTR, CDS, and 3'UTR)
#
# FPKM stands for Fragments Per Kilobase of mapped exon per Million reads.
#


import csv
import sys
import getopt
import math
import errno
import exceptions
import collections


import util

global _debug               
_debug = False                  


class Gene:
    def __init__(self):
        self.name = "name"
        self.common_name = "common_name"
        self.len = 0
        self.nfrags = 0
        self.coding = 0
        
    def __str__(self):
        rep = '{0} len: {2} \n'\
              .format(self.name, 
                      self.len )

        return rep

    pass
    # End Gene


def main(argv = None):
    """
    Prepare datafile to be used by fpkm_scatter.py.  This program
    reads from stdin and send it output to stdout.

    Input to this routine is aligned_position_stats.txt (which is the
    output of rpos_dist.py). Output from this routine is a text file
    suitable for passing to fpkm_scatter.py

    Usage:
    	calculate_fpkm.py <bowtie_out/aligned_position_stats.txt >fpkm.txt
    """
    if argv is None:
        argv = sys.argv

    mindepth = None
    try:
        opts, args = getopt.getopt(argv, "n:hd", ["help", "output="])
    except getopt.GetoptError, msg:          
        usage(msg)                         
        return 2
    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            usage()                     
            sys.exit()                  
        elif opt == '-n':
            mindepth = int(arg)                  
        elif opt == '-d':                
            global _debug               
            _debug = 1                  
        else:
            usage()                         
            return 2


    n_mapped_frags = 0
    genes = dict()
    reader = csv.DictReader(sys.stdin, delimiter='\t')
    
    for row in reader:
        gname = row["refseq"]
        if not (gname in genes):
            gene = Gene()
            gene.refseq_name = row["refseq"]
            gene.common_name = row["gene"]
            gene.len = row["genelen"]
            if row["regio"] == 'NC':
                gene.coding = 0
            else:
                gene.coding = 1
            genes[gname] = gene
        else:
            gene = genes[gname]
        gene.nfrags = gene.nfrags + 1
        n_mapped_frags = n_mapped_frags + 1

    # if the user supplied a value for mindepth (via command-line option),
    # filter out any gene that does not have at least that many
    # reads.

    if mindepth is not None:
        for key in genes.keys():
            gene = genes[key]
            avgdepth = gene.nfrags/gene.len
            if avgdepth < mindepth:
                n_mapped_frags = n_mapped_frags - gene.nfrags
                del genes[key]
                
    #
    # Calculate FPKM for each remaining gene and print the results.
    #
    hdr = "{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}"
    rec = "{0}\t{1:f}\t{2:f}\t{3:s}\t{4:d}"

    print hdr.format("tracking_id", "FPKM", "avgdepth", "gene", "coding")
    
    m = float(n_mapped_frags)/(float(1000000))
    for name in genes.keys():
        gene = genes[name]
        k = (float(gene.len)/float(1000))
        # print name, gene.len, k, m, n_mapped_frags
        fpk = (float(gene.nfrags) / k)
        fpkm = fpk / m
        avgdepth = float(gene.nfrags)/float(gene.len)

        #print name, gene.nfrags, fpkm, k, m
        print rec.format(name, fpkm, avgdepth, gene.common_name, gene.coding)
    

def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, "[-n <depth>] < <input file> "

   
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
