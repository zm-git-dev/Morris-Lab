# given a list of genes aligned against ucsc known genes, and a list
# of genes aligned agains refseq known genes, find which ones are in
# common.
#
# To do that you have to identify which ucsc genes are the same as refseq names.
#
# uc009mvx.1	0.0161775
#
# NM_009247	3.10214


import csv
import sys
import getopt
import math

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker

import numpy as np
from array import array

def main(argv = None):
    if argv is None:
        argv = sys.argv
    
    mu = 28        # default average read length
    sigma = 5      # stdev of read length 
    count = 10     # default number of reads to produce
    output = ""    # output filename (defaults to stdout)
    try:
        opts, args = getopt.getopt(argv, "hl:o:d", ["help", "output="])
    except getopt.GetoptError, msg:          
        usage(msg)                         
        return 2
    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            usage()                     
            sys.exit()                  
        elif opt == '-d':                
            global _debug               
            _debug = 1                  
        else:
            usage()                         
            return 2

    if len(args) != 3:
        usage()                         
        return 2


    ucsc_map = read_ucsc(args[0])
    refseq_map = read_refseq(args[1])
    kgXref_map = read_kgXref(args[2])


    for item in refseq_map:
        try:
            n = kgXref_map[item[0]];
            item[0] = n
        except KeyError:
            continue

    for item in refseq_map:
        print item[0], item[1]
    for item in ucsc_map:
        print item[0], item[1]
        

    return 0

def read_ucsc(filename): 
    ucsc = list()
    genemap = csv.DictReader(open(filename), delimiter='\t',
                             fieldnames=['name', 'value'])
    for row in genemap:
        try:
            ucsc.append([row['name'], float(row['value'])])
        except ValueError:
            continue

    return ucsc


def read_refseq(filename): 
    refseq = list()
    genemap = csv.DictReader(open(filename), delimiter='\t',
                             fieldnames=['name', 'value'])
    for row in genemap:
        try:
            refseq.append([row['name'], float(row['value'])])
        except ValueError:
            continue

    return refseq


def read_kgXref(filename): 

    xref = dict()
    genemap = csv.DictReader(open(filename), delimiter='\t')

    for row in genemap:
        try:
            xref[row['mRNA']] = row['#kgID']
        except ValueError:
            continue

    return xref






if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
