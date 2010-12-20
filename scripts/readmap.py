# read the refseq map into a dictionary

# this is a poor-man's in-memory database
#1	DV038825	243377
#
#

import csv
import sys
import getopt
import math

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker

import numpy as np
from array import array



# main() takes an optional 'argv' argument, which allows us to call it
# from the interactive Python promp.

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

    platform_map = read_map(args[0])
    
    ref_map = read_unigene2refseq(args[1])

    process_sample_file(args[2], platform_map, ref_map)

    return 0


def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <refseq.map> "




# Read the platform map
# 1	DV038825	243377
# 3	BM877548	76748
# 6	BC057190	104625
#
def read_map(reference): 
    print "read_map"
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


# Read the refseq_reflink.txt
#name	product	mrnaAcc	protAcc	geneName	prodName	locusLinkId	omimId
# MIR375		NR_035777		390601	0	100316383	0
# Os08g0544800	hypothetical protein	NM_001068944	NP_001062409	235443	97	4346214	0
# Os04g0476700	hypothetical protein	NM_001059617	NP_001053082	265619	97	4336157	0
# Os01g0170300	hypothetical protein	NM_001048674	NP_001042139	233119	97	4327633	0
#
def read_unigene2refseq(reflink_file):
    print "read_unigene2refseq"
    refmap = dict()
    refmap_reader = csv.DictReader(open(reflink_file), delimiter='\t')
    for row in refmap_reader:
        refmap[row['locusLinkId']] = row['mrnaAcc']

    return refmap


# 1 	0.1771207
# 2 	-0.4996132
# 3 	0.07726434
# 4 	-0.3684799
def process_sample_file(sample_file, platform_map, ref_map):

    refmap = dict()
    sample_reader = csv.DictReader(open(sample_file), delimiter='\t', fieldnames=['#id', '#abundance'])

    for row in sample_reader:
        try:
            entrez = platform_map[int(row['#id'])]
            gene = ref_map[entrez]
            abundance = math.pow(2, float(row['#abundance']))
            print gene, abundance
            continue
        except ValueError:
            continue
        except KeyError:
            continue
                


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
