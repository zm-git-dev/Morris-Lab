# Find where  a read is mapped w.r.t. the coding start region.

# Works with a wide format knowngene file, NOT GTF file

# Start with the list of known genes.  Read the whole thing into memory - it's not that big.



import csv
import sys
import getopt
import math
import pprint

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker

import numpy as np
from array import array
from itertools import ifilter

import gnote
import genexref
import util

from Bio.SeqFeature import FeatureLocation

from Bio import SeqIO
from Bio import SeqFeature
from BCBio import GFF
from BCBio.GFF import (GFF3Writer, GFFExaminer, GFFParser, DiscoGFFParser)

import pysam



class Gene:
    pass

class Exon:
    pass


# main() takes an optional 'argv' argument, which allows us to call it
# from the interactive Python promp.

def main(argv = None):
    if argv is None:
        argv = sys.argv
    
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

    if len(args) != 2:
        usage()                         
        return 2

    known_genes = args[0]
    accepted_hits = args[1]

    chromosomes = dict()


    #bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	id	name2	cdsStartStat	cdsEndStat	exonFrames
    #608	NR_003519	chr17	-	3064317	3084183	3084183	3084183	8	3064317,3076499,3077292,3077687,3078999,3079354,3080868,3083967,	3064460,3076782,3077446,3077834,3079138,3079591,3081034,3084183,	172717	Pisd-ps2	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,
    #9	NM_134123	chr17	+	3114971	3198859	3115388	3198210	20	3114971,3145074,3148535,3159162,3162958,3164146,3167962,3169326,3171102,3177092,3177610,3178154,3185063,3185760,3187588,3190146,3191816,3193013,3195777,3196759,	3115418,3145158,3148580,3159324,3163112,3164277,3168139,3169406,3171220,3177224,3177723,3178348,3185164,3185874,3187745,3190280,3191961,3193082,3195993,3198859,	277518	Scaf8	cmpl	cmpl	0,0,0,0,0,1,0,0,2,0,0,2,1,0,0,1,0,1,1,1,


    # read in the known genes.
    #
    try:
        in_handle = None
        in_handle = open(known_genes)

        previous = Gene()
        previous.end = 0
        g_reader = csv.DictReader(in_handle, delimiter='\t')
        for g_row in g_reader:
            # parse exons.
            es_reader = csv.reader([g_row["exonStarts"]], delimiter=',')
            ee_reader = csv.reader([g_row["exonEnds"]], delimiter=',')
            estart_list = es_reader.next()
            eend_list = ee_reader.next()
            exon_list = list()
            for i in range(int(g_row["exonCount"])):
                #print i, int(estart_list[i]), int(eend_list[i])
                exon_list.append([int(estart_list[i]), int(eend_list[i])])

            # parse .
            es_reader = csv.reader([g_row["exonStarts"]], delimiter=',')
            ee_reader = csv.reader([g_row["exonEnds"]], delimiter=',')
            ec_reader = csv.reader([g_row["exonFrames"]], delimiter=',')
            estart_list = es_reader.next()
            eend_list = ee_reader.next()
            ecnd_list = ec_reader.next()
            exon_list = list()
            for i in range(int(g_row["exonCount"])):
                #print i, int(estart_list[i]), int(eend_list[i])
                exon = Exon()
                exon.start = int(estart_list[i])
                exon.end = int(eend_list[i])
                exon.frame = int(ecnd_list[i])
                exon_list.append(exon)

            gene = Gene()
            gene.start = int(g_row["txStart"])
            gene.end = int(g_row["txEnd"])
            gene.exons = exon_list
            gene.cdsStart = int(g_row["cdsStart"])
            gene.cdsEnd = int(g_row["cdsEnd"])
            gene.name = g_row["name"]
            if  g_row["strand"] == "-":
                gene.strand = -1
            else:
                gene.strand = +1
                
            assert gene.end >= previous.end
            previous = gene
            

            chrom = g_row["chrom"]
            if chrom not in chromosomes:
                chromosomes[chrom]= list()
            chromosomes[chrom].append(gene)

            
    finally:
        if in_handle is not None:
            in_handle.close()






    # process the reads
    #
    try:
        samfile = None
        samfile = pysam.Samfile( accepted_hits, "r" )

        unmatched = 0
        unmatched_gene = 0
        unmatched_strand = 0
        unmatched_exon = 0
        nreads = 0
        warnings = dict()
        for alignedread in samfile.fetch():
            nreads += 1
            #print "read aligned read from SAM file"
            matched = False
            matched_gene = False
            matched_strand = False
            matched_exon = False
            chrom = samfile.getrname(alignedread.rname)
            if chrom not in chromosomes:
                if chrom not in warnings:
                    print 'Warning:   Could not find known genes for chromosome "%s"' % (chrom)
                    warnings[chrom] = True
                continue
            previous = Gene()
            previous.end = 0
            for gene in chromosomes[chrom]:

                assert gene.end >= previous.end
                previous = gene

                #
                # Assume the known gene list is sorted by ending
                # position.  If the current read starts at a point
                # past the end of the current gene, it will also be
                # past every subsequent gene.  Just skip all the
                # follwing genes because there is no way they can
                # match
                #
                if alignedread.pos > gene.end:
                    break;
                
                #
                # check whether this read overlaps this gene...
                #
                if overlap(alignedread, gene):
                    #print "read %s overlaps gene %s" % (alignedread.qname, gene.name)
                    matched_gene = True
                    if strands_match(alignedread, gene):
                        #print "\tstrands match"
                        matched_strand = True
                        exonlist = gene.exons
                        for exon in exonlist:
                            # check whether this read overlaps this exon...
                            if overlap(alignedread, exon):
                                #print "\t\tread overlaps exon"
                                if matched == True:
                                    print "READ MATCHED MORE THAN ONE EXON!"
                                    print "[%d, %d] (%d) %d " % (alignedread.pos, alignedread.pos+alignedread.rlen, alignedread.rlen, alignedread.is_reverse )
                                    print "[%d, %d] (%d) " % (exon.start, exon.end, exon.end-exon.start)
                                    print "\n"
                                matched = True
                                matched_exon = True
                                break

                # if the read matched this gene, stop looking.
                if matched:
                    break

            if not matched:
                unmatched += 1
                if not matched_gene:
                    unmatched_gene += 1
                elif not matched_strand:
                    unmatched_strand += 1
                elif not matched_exon:
                    unmatched_exon += 1
            

    finally:
        if samfile is not None:
            samfile.close()

    print "\n%d reads processed.\n" % (nreads)
    print "%d matched an exon within a gene.\n" % (nreads - unmatched)
    print "%d did not match." % (unmatched)
    print "\t%d did not match any gene" % (unmatched_gene)
    print "\t%d were on the wrong strand" % (unmatched_strand)
    print "\t%d did not match an exon" % (unmatched_exon)

    return 0


# determine whether the exon "contains" the read
#
def overlap(read, feature):
    assert feature.start < feature.end 
    return (read.pos >=  feature.start and read.pos+read.rlen <= feature.end)
                
def strands_match(read, gene):
    assert gene.strand == -1 or gene.strand == 1
    return (read.is_reverse and gene.strand == -1) or (not read.is_reverse and gene.strand == 1)
                

def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <reference_tbl.csv> <alignments.bam> "

   
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
