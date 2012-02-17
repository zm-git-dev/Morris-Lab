#! /usr/bin/python
# compare_depth.py exp1/aligned_position_stats.txt exp2/aligned_position_stats.txt >output
#
# This script takes two aligned_position_stats.txt files from two
# different experiments and prepares a table that compares the read
# depth for the union of all genes covered by the two experiments.
#
#
# The input file for this script typically comes from rpos_dist.py
#
# input format looks like this:
#    read name       refseq              gene     genelen    region    %    gene%
# ILLUMINA-D43DB5  NM_001164655    9530053A07Rik  3083        CDS    39.80  39.16
#
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
        self.depth = 0
        self.coding = 0
        
    def __str__(self):
        rep = '{0} len: {2} \n'\
              .format(self.name, 
                      self.len )
        return rep

    # End Gene

class GeneFile:
    def __init__(self):
        self.name = ""
        self.genes = dict()
        self.mapped_frags = 0
        
    def __str__(self):
        rep = '{0} len: {2} \n'\
              .format(self.name, 
                      self.len )

        return rep

    def read_genefile(self, fname):
        self.name = fname
        in_handle = None
        try:
            in_handle = open(fname)
            reader = csv.DictReader(in_handle, delimiter='\t')
            for row in reader:
                gname = row["refseq"]
                if not (gname in self.genes):
                    gene = Gene()
                    gene.refseq_name = row["refseq"]
                    gene.common_name = row["gene"]
                    gene.len = row["genelen"]
                    if row["regio"] == 'NC':
                        gene.coding = 0
                    else:
                        gene.coding = 1
                    self.genes[gname] = gene
                else:
                    gene = self.genes[gname]
                gene.depth = gene.depth + 1
                self.mapped_frags = self.mapped_frags + 1

        finally:
            if in_handle is not None:
                in_handle.close()

    # End GeneFile


def main(argv = None):
    if argv is None:
        argv = sys.argv

    mindepth = None
    genename = None
    try:
        opts, args = getopt.getopt(argv, "g:n:hd", ["help", "output="])
    except getopt.GetoptError, msg:          
        usage(msg)                         
        return 2
    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            usage()                     
            sys.exit()                  
        elif opt == '-n':
            mindepth = int(arg)                  
        elif opt == '-g':
            genename = arg                  
        elif opt == '-d':                
            global _debug               
            _debug = 1                  
        else:
            usage()                         
            return 2

    genef1 = GeneFile()
    genef2 = GeneFile()

    genef1.read_genefile(args[0])
    genef2.read_genefile(args[1])

    hdr = "{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\t{5:s}"
    rec = "{0:s}\t{1:s}\t{2:d}\t{3:d}\t{4:d}\t{5:d}"

    print hdr.format("NM_name", "common_name", "gene_len", "mfrag1", "mfrag2", "coding")

    if genename != None:
        combined = set([genename])
    else:
        geneset1 = set(genef1.genes.keys())
        geneset2 = set(genef2.genes.keys())
        combined = geneset1.union(geneset2)
    for k in combined:
        depth1 = depth2 = 0

        if k in genef1.genes:
          gene = genef1.genes[k]
          depth1 = genef1.genes[k].depth

        if k in genef2.genes:
          gene = genef2.genes[k]
          depth2 = genef2.genes[k].depth

        #print name, gene.depth, fpkm, k, m
        print rec.format(k, gene.common_name, int(gene.len), depth1, depth2, gene.coding)


def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, "[-n <depth>] < <input file> "

   
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
