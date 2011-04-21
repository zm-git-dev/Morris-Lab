# Morris Lab Tools
#
# Some usefule tools (mostly file parsing utilities) that are used by
# a number of other scripts.

import csv
import sys
import math
import pprint

import numpy as np
from array import array

import genexref
import util

from Bio.SeqFeature import FeatureLocation

from Bio import SeqIO
from Bio import SeqFeature

import pysam

class Gene:
    pass

class CDS:
    pass

class Exon:
    def len(self):
       return self.end - self.start 
    pass


class FileFormatException(Exception):
    def __init__(self, file, field):
        self.filename = file
        self.fieldname = field
    def __str__(self):
         return repr(self.fieldname)
     

# Read a .TXT file withone entry for each known gene.  The file format
# is TXT and the columns are arranged in the following order:
#
# #bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount
# #exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames
#
# http://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Mouse&db=mm9&hgta_group=genes&hgta_track=refGene&hgta_table=0&hgta_regionType=genome
#
                
def read_knowngenes(known_genes, genename = None):
    chromosomes = dict()
    genecount = 0
    badframecount = 0


    # Files parsed by this routine must have the following columns.
    # The order of the columns does not matter.  I dont try to
    # validate the data more thoroughly at this time.  If you got this
    # file from the right table at UCSC genome browser you should be
    # OK.
    #
    HeaderFields = [ "name", "chrom", "strand", "txStart", "txEnd",
                     "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                     "id", "name2", "cdsStartStat", "cdsEndStat", "exonFrames" ]

    # read in the known genes.
    #
    try:
        in_handle = open(known_genes)
        
        g_reader = csv.DictReader(in_handle, delimiter='\t')

        for field in HeaderFields:
            if not field in g_reader.fieldnames:
                raise FileFormatException(known_genes, field)

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

    except FileFormatException as (e):
        print "file format exception: file %s must have a column called \"%s\"" % (e.filename, e.fieldname)
        
    finally:
        if in_handle is not None:
            in_handle.close()

    return chromosomes
    # END read_knownegens()

    

