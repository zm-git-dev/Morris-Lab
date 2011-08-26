import csv
import sys
import getopt
import math
import pprint

import numpy as np
from array import array
from itertools import ifilter

#import gnote
import genexref
import util

from Bio.SeqFeature import FeatureLocation

from Bio import SeqIO
from Bio import SeqFeature
#from BCBio import GFF
#from BCBio.GFF import (GFF3Writer, GFFExaminer, GFFParser, DiscoGFFParser)

import pysam
import mltools

global _debug               
_debug = False                  

# main() takes an optional 'argv' argument, which allows us to call it
# from the interactive Python promp.

def main(argv = None):
    if argv is None:
        argv = sys.argv
    
    gene_name = None
    match_limit = None
    range_str = None
    print_genelist = False
    
    try:
        opts, args = getopt.getopt(argv, "hdLl:g:r:", ["help", "output="])
    except getopt.GetoptError, msg:          
        usage(msg)                         
        return 2
    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            usage()                     
            sys.exit()                  
        elif opt == '-d':                
            global _debug               
            _debug = True                  
        elif opt == '-L':  # print a gene list and exit.
            print_genelist = True
        elif opt == '-r':  # range of nucleotide to look at.
                           # takes a range as "n:m" where
                           # n < m
                           # -r -200:200
            range_str = arg
        elif opt == '-g':  # limit the match to one particular gene
                           # takes a gene name as argument;
                           # -g NM_29391293
            gene_name = arg
        elif opt == '-l':  # limit the number of matching reads
                           # takes an integrer as argument
                           # stop processing after matching n reads to exons.
            match_limit = int(arg)
        else:
            usage()                         
            return 2

    if len(args) != 2:
        usage()                         
        return 2

    gene_file = args[0]
    reads_file = args[1]

    knowngenes = mltools.read_knowngenes(gene_file, gene_name)
    collapse_genome(knowngenes)
    # export_transcriptome(knowngenes)
    
    return 0


def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <reference_tbl.csv> <alignments.bam> "


# Output in FASTA format suitable for use by IGV
# First collapse the genome then output the FASTA data.
#
# By "collapse genome" I mean eliminate introns.  For a given gene,
# leave the first exon in place and shift every subsequent exon toward
# the first exon to eliminate gaps.  the same algorith applies to
# genes on the reverse strand, only those will be shifted in the
# opposit direction.
#
def collapse_genome(knowngenes):
    for k, v in chromosomes:
        for genes in v:
            if (gene.strand == 1):
                # forward strand
                for i in range(2, len(gene.exons)):
                    exon = gene.exons[i] 
                    previous = gene.exons[i-1]
                    intron_len = exon.start - previous.end
                    exon.start -= intron_len
                    exon.end -= intron_len

            else:
                # reverse strand
                for i in range(len(gene.exons)-1, 1):
                    exon = gene.exons[i]
                    previous = gene.exons[i+1]
                    intron_len = previous.start - exon.end
                    exon.start += intron_len
                    exon.end += intron_len


   
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
