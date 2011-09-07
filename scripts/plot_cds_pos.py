# plot CDS position
# uses information from the rpos_dist.py
# plots mod%3 of CDS position.
# python /home/csw/Morris-Lab/scripts/test.py  /usr/local/share/genome/mm9/mm9_refseq_knowngenes.txt tophat_out/overlaps_exons_uniq.bed | awk '{print $4}'
import csv
import re
import sys
import getopt
import math
import pprint

import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#import matplotlib.ticker as ticker

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

global _debug               
_debug = False                  

import fileinput

def graph_startpos(startpos, range_str):
    """
    Make a histogram of normally distributed random numbers and plot the
    analytic PDF over it
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if range_str != None:
        parts = range_str.split(':')
        low = int(parts[0])
        high = int(parts[1])
        n, bins, patches = plt.hist(startpos, high-low, range=(low, high), histtype='bar')
    else:
        n, bins, patches = plt.hist(startpos, 900, histtype='bar')

        
    # the histogram of the data with histtype='step'
    #n, bins, patches = plt.hist(startpos, 50, normed=1, histtype='stepfilled')
    plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

    ax.set_xlabel('CDS Position [nt from start]')
    ax.set_ylabel('Read 5\' ends')
    ax.set_title(r'$\mathrm{Count\ of\ read\ alignments\ per\ gene}$')

    plt.show()



def main(argv = None):
    if argv is None:
        argv = sys.argv
    
    gene_name = None
    match_limit = None
    range_str = None
    print_genelist = False
    
    try:
        opts, args = getopt.getopt(argv, "hdg:r:", ["help", "output="])
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
        elif opt == '-r':  # range of nucleotide to look at.
                           # takes a range as "n:m" where
                           # n < m
                           # -r -200:200
            range_str = arg
        elif opt == '-g':  # limit the match to one particular gene
                           # takes a gene name as argument;
                           # -g NM_29391293
            gene_name = arg
        else:
            usage()                         
            return 2


    poslist = list()
    for line in sys.stdin:
        if not re.match("^#", line):
            [rname, nmname, cname, cds_pos, rest] = line.split(None, 4)
            poslist.append(int(cds_pos))
            print cds_pos
        
    graph_startpos(poslist, range_str)

    return 0



def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <reference_tbl.csv> <alignments.bam> "

   
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
