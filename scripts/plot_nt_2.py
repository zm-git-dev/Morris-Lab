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


class Gene:
    pass

class CDS:
    pass

class Exon:
    def len(self):
       return self.end - self.start 

    pass

class Stats:
    pass

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

    knowngenes = read_knowngenes(gene_file, gene_name)
    poslist = process_reads(reads_file, knowngenes, match_limit)

    if print_genelist:
        print "# common_name\trefseq_name\tread_count"
        for chrom in knowngenes:
            for gene in knowngenes[chrom]:
                print "%s\t%s\t%d" % (gene.common_name, gene.name, gene.readcount)
    else:
        graph_startpos(poslist, range_str)

    return 0


# determine whether the range "contains" the read
#
def overlap(read, feature):
    assert feature.start < feature.end
    assert feature.start > 0
    assert read.pos > 0 
    return (read.pos >=  feature.start and read.pos+read.rlen <= feature.end)
                
def strands_match(read, gene):
    assert gene.strand == -1 or gene.strand == 1
    return (read.is_reverse and gene.strand == -1) or (not read.is_reverse and gene.strand == 1)
                
def read_knowngenes(known_genes, genename = None):
    chromosomes = dict()
    genecount = 0
    badframecount = 0

    # read in the known genes.
    #
    try:
        in_handle = open(known_genes)
        
        g_reader = csv.DictReader(in_handle, delimiter='\t')
        for g_row in g_reader:
            if genename != None and genename != g_row["name"]:
                continue
            gene = Gene()
            gene.start = int(g_row["txStart"])
            gene.end = int(g_row["txEnd"])
            gene.cds = CDS()
            gene.cds.start = int(g_row["cdsStart"])
            gene.cds.end = int(g_row["cdsEnd"])
            gene.readcount = 0
            gene.name = g_row["name"]
            gene.common_name = g_row["name2"]

            if gene.cds.start == gene.cds.end:
                # it is a non-coding RNA.
                # don't add to the gene list.
                continue
            
            gene.name = g_row["name"]
            if  g_row["strand"] == "-":
                gene.strand = -1
            else:
                gene.strand = +1

            
            es_reader = csv.reader([g_row["exonStarts"]], delimiter=',')
            ee_reader = csv.reader([g_row["exonEnds"]], delimiter=',')
            ef_reader = csv.reader([g_row["exonFrames"]], delimiter=',')
            estart_list = es_reader.next()
            eend_list = ee_reader.next()
            efrm_list = ef_reader.next()
            exon_list = list()
            for i in range(int(g_row["exonCount"])):
                exon = Exon()
                exon.start = int(estart_list[i])
                exon.end = int(eend_list[i])
                exon.frame = int(efrm_list[i])

                exon_list.append(exon)


            # calculate how far from the start of the coding region each exon begins.
            
            cds_i = 0
            cds_e = 0
            if (gene.strand == 1):
                # forward strand

                #Identify which exon contains the CDS.
                for i in range(0,len(exon_list)):
                    exon = exon_list[i]
                    if (gene.cds.start >= exon.start and gene.cds.start <= exon.end):
                        cds_i = i
                    if (gene.cds.end >= exon.start and gene.cds.end <= exon.end):
                        cds_e = i
                        
                                        
                exon = exon_list[cds_i]
                exon.pos_cds = exon.start - gene.cds.start
                for i in range(cds_i+1,len(exon_list)):
                    exon = exon_list[i]
                    exon.pos_cds = exon_list[i-1].pos_cds + exon_list[i-1].len()

                for i in range(cds_i-1,-1,-1):
                    exon = exon_list[i]
                    exon.pos_cds = exon_list[i+1].pos_cds - exon.len()
            else:
                # reverse strand

                #Identify which exon contains the CDS.
                for i in range(0,len(exon_list)):
                    exon = exon_list[i]
                    if (gene.cds.end >= exon.start and gene.cds.end <= exon.end):
                        cds_i = i
                    if (gene.cds.start >= exon.start and gene.cds.start <= exon.end):
                        cds_e = i

                exon = exon_list[cds_i]
                exon.pos_cds = gene.cds.end - exon.end
                for i in range(cds_i+1, len(exon_list)):
                    exon = exon_list[i]
                    exon.pos_cds = exon_list[i-1].pos_cds - exon.len()

                for i in range(cds_i-1,-1,-1):
                    exon = exon_list[i]
                    exon.pos_cds = exon_list[i+1].pos_cds + exon_list[i+1].len()

#                 print "%s %d" % (gene.name, gene.strand)
#                 for i in range(0,len(exon_list)):
#                     exon = exon_list[i]
#                     if (i == cds_i):
#                         print "*%d\t%d\t%d\t%d (%d) (expect %d)" % (exon.start, exon.end, exon.len(), exon.pos_cds, exon.pos_cds%3, exon.frame)
#                     else:
#                         print " %d\t%d\t%d\t%d (%d) (expect %d)" % (exon.start, exon.end, exon.len(),  exon.pos_cds, exon.pos_cds%3, exon.frame)
#                 print "\n\n"
                
            gene.exons = exon_list
                
            chrom = g_row["chrom"]
            if chrom not in chromosomes:
                chromosomes[chrom]= list()
            chromosomes[chrom].append(gene)
            genecount += 1

            
    finally:
        if in_handle is not None:
            in_handle.close()

    return chromosomes
    # END read_knownegens()

    

# process the reads
#
def process_reads(reads_file, knowngenes, match_limit):

    try:
        poslist = list()
        samfile = None
        samfile = pysam.Samfile( reads_file, "r" )

        unmatched = 0
        unmatched_gene = 0
        unmatched_strand = 0
        unmatched_exon = 0
        nmatched = 0
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
            if chrom not in knowngenes:
                if chrom not in warnings:
                    print >> sys.stderr, 'Warning:   Could not find known genes for chromosome "%s"' % (chrom)
                    warnings[chrom] = True
                continue
            for gene in knowngenes[chrom]:

                #
                # check whether this read overlaps this gene...
                #
                if overlap(alignedread, gene):
                    matched_gene = True
                    if strands_match(alignedread, gene):
                        matched_strand = True
                        exonlist = gene.exons
                        for exon in exonlist:
                            # check whether this read overlaps this exon...
                            if overlap(alignedread, exon):
                                if matched == True:
                                    print >> sys.stderr, "READ MATCHED MORE THAN ONE EXON!"
                                    print >> sys.stderr, "%s\t[%d, %d] (%d) %d " % (chrom, alignedread.pos, alignedread.pos+alignedread.rlen, alignedread.rlen, alignedread.is_reverse )
                                    print >> sys.stderr, "\t[%d, %d] (%d) " % (exon.start, exon.end, exon.end-exon.start)
                                    print >> sys.stderr, "\n"
                                gene.readcount += 1
                                matched = True
                                matched_exon = True
                                if gene.strand == -1:
                                    pos = exon.pos_cds + (exon.end - alignedread.pos)
                                else:
                                    pos = exon.pos_cds + (alignedread.pos - exon.start)
                                poslist.append(pos)

                                global _debug
                                if _debug:
                                    print "%s matched %s:%s" % (alignedread.qname, gene.name, gene.common_name)
                                
                                break # for exon in exonlist

                # if the read matched this gene, stop looking.
                if matched:
                    break   #  for gene in knowngenes[chrom]
                    
            if not matched:
                unmatched += 1
                if not matched_gene:
                    unmatched_gene += 1
                elif not matched_strand:
                    unmatched_strand += 1
                elif not matched_exon:
                    unmatched_exon += 1
            else:
                nmatched += 1

            # temporarily limit number of reads to get a sense of whather this is working.
            
            if match_limit != None and nmatched >= match_limit:
                break

    finally:
        if samfile is not None:
            samfile.close()

#     print "\n%d reads processed.\n" % (nreads)
#     print "%d matched an exon within a gene.\n" % (nreads - unmatched)
#     print "%d did not match." % (unmatched)
#     print "\t%d did not match any gene" % (unmatched_gene)
#     print "\t%d were on the wrong strand" % (unmatched_strand)
#     print "\t%d did not match an exon" % (unmatched_exon)

    return poslist
    # END process_reads()


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








def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <reference_tbl.csv> <alignments.bam> "

   
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
