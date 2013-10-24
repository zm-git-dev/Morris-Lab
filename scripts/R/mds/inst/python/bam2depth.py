#!/usr/bin/python -O
##
## Annotate
## Author: Chris Warth
##
## bam2depth - calculate the depth of reads across gene annotations as
## specificed in a GTF file.
##
## If you have a BAM file and a gene model in a GTF file (like RefGenes or the knownGenes model from UCSC)
## this routine will output the depth of reads at each nucleotide of a 

##

from __future__ import print_function


import sys
import getopt
import csv
import collections
import subprocess
import os

""" Print a usage message
"""
def usage(msg = None):
    if msg is not None:
        print(msg)
    print(__name__, " <gtf_file> <bam files>")

def run_command(command):
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')

# list flattening code borrowed from:
# http://stackoverflow.com/a/406822/1135316
def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def processGene(bamfiles, transID, chrom, strand, exons):
    # if on the reverse strand, exons are listed in reverse order
    if strand == '-':
        exons = list(reversed(exons))

    # keep track of the transcript position of each exon.
    exonPos = 0
    exonNum = 1
    for exon in exons:
        cmd = ["samtools", "depth", bamfiles, "-r", "%s:%d-%d" % (chrom, exon[0], exon[1])]
        cmd = flatten(cmd)
        
        # print(cmd, file=sys.stderr)
        for line in run_command(cmd):
            line = line.rstrip("\n")
            # parse the line to get the genomic position so we can calculate the transcript position
            depth =line.split()
            genomePos = int(depth[1])
            if strand == '+':
                transPos = exonPos + genomePos - exon[0] + 1
            else:
                transPos = exonPos + exon[1] - genomePos + 1
            print("%s\t%d\t%d\t%s" % (transID, transPos, exonNum, line)),

        exonPos += abs(exon[1]-exon[0])+1
        exonNum += 1


def processGTF(gtffile, bamfiles, geneList):
    """
    return a dictionary to sequences.  Use the sequence identifier as the key.
    """

    ## damn, there is no way to debug a python
    ## program that expects file input from stdin.  When someone asks about this
    ## egregious behavior the replies are all the same - change your
    ## program to read from a file.  I would rather have this program
    ## take the gtf file on stdin, but if I do that python cannot
    ## debug it.

    hdr = "symbol\ttransPos\texonNumber\tchrm\tgenomePos"
    for bfile in bamfiles:
        tmp = os.path.splitext(os.path.basename(bfile))[0]
        hdr += "\t"+tmp
    print(hdr)
    
    transID = None
    exons = []
    with open(gtffile,'rb') as gtfin:
        lineno = 0
        gtfin = csv.reader(gtfin, delimiter='\t')
        for row in gtfin:
            lineno += 1
            # each row of the GTF should have 9 fields.  Probably should relax this
            # a little to allow more than 9, but not fewer.
            assert len(row) == 9, "line %d: number of fields (%r) is not 9." % (lineno, len(row))
            
            # samtools depth will barf if it encounters a chromosome name with an underscore
            # so filter those out.
            #
            # $ samtools depth test.bam -r chr17_c:195721-200597
            # [bam_parse_region] fail to determine the sequence name.

            if row[0].find('_') != -1:
                continue
            
            # The last field is a series of name-value pairs, separated by ';'
            annotations = row[8].split(";")
            attr = {}
            for i in annotations:
                (name, value) = i.split()
                attr[name] = value.strip('"')

            if (not geneList or (attr["transcript_id"] in geneList)):
                # entries with identical transcript IDs come from the same gene.
                if (transID != attr["transcript_id"]):
                    if exons:
                        processGene(bamfiles, transID, chrom, strand, exons)
                    transID = attr["transcript_id"]
                    chrom = row[0]
                    strand = row[6]
                    exons = []

                start = int(row[3])
                end = int(row[4])
                exons.append((start, end))

                # assert that the strand and chromosome do not change halfway through the gene annotation
                assert strand == row[6], "line %d, Gene %s: inconsistent orientation!\n%s %s" % ( lineno, transID, strand, row[6])
                assert chrom == row[0], "line %d, Gene %s: inconsistent chromosome!\n%s %s" % ( lineno, transID, chrom, row[0])


        # After the GTF file is exhausted, output any remaining data from a reverse-strand gene.
        if exons:
            processGene(bamfiles, transID, chrom, strand, exons)
            transID = attr["transcript_id"]
                    



def main(argv = None):
    genelist = []
    if argv is None:
        argv = sys.argv
    try:
        opts, args = getopt.getopt(argv, "hg:", ["help", "output="])
    except getopt.GetoptError, msg:          
        usage(msg)                         
        return 2
    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            usage()                     
            sys.exit()                  
        elif opt in ("-g", "--gene"):      
            genelist.append(*arg.split(","))
        elif opt == '--':  # end of options
            break
        else:
            usage()                         
            return 2

    gtfFile = args[0]
    bamFiles = args[1:]

    processGTF(gtfFile, bamFiles, genelist)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

        

