# make a plot of ribosome position vs. frequency.
# This is one of the graphs we discussed at the hutch.
#
#
# To visualize the results from this routine,
#
# $ python  /home/csw/Morris-Lab/scripts/test.py  /home/csw/Morris-Lab/ref-seq/mm9/refseq_knowngenes.txt overlaps_exons_uniq.bed >out.txt
#
# $ sed '/^#/d' <out.txt | sort -k 6n | uniq -c -f 5 | awk '{print $1, $7}' >out2.txt
#
# $ gnuplot <<EOF
# set xlabel "Position on gene"
# set ylabel "Number of reads"
# set title "Read position vs. Read count"
# plot  "out2.txt" using 2:1 title '' with lines
# EOF
#

import csv
import sys
import getopt
import math
import pprint
import errno
import exceptions
import collections


#import gnote
import genexref
import util

global _debug               
_debug = False                  

class CDS:
    def __init__(self):
        self.start = 0
        self.end = 0
    def __str__(self):
        rep = '{0} {1}'.format(self.start, self.end)
        return rep
    pass


class Gene:
    def __init__(self):
        self.exons = list()
        self.start = 0
        self.end = 0
        self.cds = CDS()
        self.readcount = 0
        self.name = "name"
        self.common_name = "common_name"
        self.len = 0


    def __str__(self):
        rep = '{0} {1!s} len: {2} ({3} {4} {5})  CDS: {6!s}\n'\
              .format(self.name, self.strand,
                      self.len, self.len_utr_5(), self.len_cds(), self.len_utr_3(),
                      self.cds)
        rep += ' {0:^10} {1:^10} {2:^5} {3:^9s} {4:^5s} {5:^5s} {6:^5s}\n'\
               .format("start", "end", "len", "cds", "sog", "<cds", "cds>");
        icds_start = self.exon_containing(self.cds.start)
        icds_end = self.exon_containing(self.cds.end)

        # If the gene is on the reverse strand, swap the CDS 
        if (self.strand == -1):
            icds_start, icds_end = icds_end, icds_start

        for i in xrange(0,len(self.exons)):
            exon = self.exons[i]
            if (i == icds_start):
                star = ">"
            elif (i == icds_end):
                star = "<"
            else:
                star = ""
            rep += "{0:1}{1:s}\n".format(star, exon)

        return rep


    # Identify which exon contains the CDS.  Return -1 if there is no
    # coding region (noncoding rna)
    def exon_containing(self, pos):
        cds_i = -1
        for i in xrange(0,len(self.exons)):
            exon = self.exons[i]
            if (self.strand == 1):
                if (pos >= exon.start and pos <= exon.end):
                    cds_i = i
            else:
                if (pos >= exon.start and pos <= exon.end):
                    cds_i = i
        return cds_i

    #
    # Calculate how far the start of each exon is from the ends of the coding region.
    #
    def annotate_cds(self):
        self.annotate_cds_start()
        self.annotate_cds_end()

    # Calculate how far the start of each exon is from the beginnging
    # of the coding region.
    #
    def annotate_cds_start(self):
        # find the exon that contains the position
        
        if (self.strand == 1):
            # forward strand
            cds_i = self.exon_containing(self.cds.start)
            exon = self.exons[cds_i]
            exon.cds.start = exon.start - self.cds.start
        else:
            # reverse strand
            cds_i = self.exon_containing(self.cds.end)
            exon = self.exons[cds_i]
            exon.cds.start = self.cds.end - exon.end

        # proceed in the forward direction, annotating all
        # the exons *after* the CDS
        #
        for i in xrange(cds_i+1,len(self.exons)):
            self.exons[i].cds.start = self.exons[i-1].cds.start + self.exons[i-1].len()
 
        # proceed in the reverse direction, annotating 
        # the exons *before* the CDS
        #
        for i in xrange(cds_i-1,-1,-1):
            self.exons[i].cds.start = self.exons[i+1].cds.start - self.exons[i].len()


    # Calculate how far the start of each exon is from the end
    # of the coding region.
    #
    def annotate_cds_end(self):
        if (self.strand == 1):
            # forward strand
            cds_i = self.exon_containing(self.cds.end)
            exon = self.exons[cds_i]
            exon.cds.end = exon.start - self.cds.end
        else:
            # reverse strand
            cds_i = self.exon_containing(self.cds.start)
            exon = self.exons[cds_i]
            exon.cds.end =  self.cds.start- exon.end
        
        # proceed in the forward direction, annotating all
        # the exons *after* the CDS
        #
        for i in xrange(cds_i+1,len(self.exons)):
            self.exons[i].cds.end = self.exons[i-1].cds.end + self.exons[i-1].len()
 
        # proceed in the reverse direction, annotating 
        # the exons *before* the CDS
        #
        for i in xrange(cds_i-1,-1,-1):
            self.exons[i].cds.end = self.exons[i+1].cds.end - self.exons[i].len()


    def len_cds(self):
        return self.len - (self.len_utr_5() + self.len_utr_3())

    def len_utr_5(self):
        return abs(self.exons[0].cds.start)

    def len_utr_3(self):
        nexon = len(self.exons)-1
        return self.exons[nexon].cds.end + self.exons[nexon].len()

    pass
    # End Gene


class Exon:
    def __init__(self):
        self.start = 0
        self.end = 0
        self.frame = +1
        self.cds = CDS()

    def len(self):
       return abs(self.end - self.start)

    def __str__(self):
        if (self.cds.start%3 != self.frame):
            oof = "*"
        else:
            oof = ""

        rep = "{0:10d} {1:10d} {2:5d} {3:5d} ({4}{5:1}) {6:5d} {7:5d} {8:5d}"

        return rep.format(self.start, self.end, self.len(), self.cds.start, \
                          self.cds.start%3, oof, self.pos_gene, self.cds.start, self.cds.end) 
    pass
    # End Exon

class Stats:
    def __init__(self):
        self.reads_processed = 0
        self.reads_mapped = 0
        self.genes_covered = dict()
        self.multi_locus_read = 0
        self.region_5prime = 0
        self.region_3prime = 0
        self.region_coding = 0
        self.region_noncoding = 0

    def __str__(self):
        rep = ""
        rep += "reads processed:   {0:d}\n".format(self.reads_processed)
        rep += "reads mapped:      {0:d}\n".format(self.reads_mapped)
        rep += "reads mapped to \n"
        rep += "multi_locus genes: {0:d}\n".format(self.multi_locus_read)

        nreads = 0
        for v in self.genes_covered.values():
            nreads += v
        ngenes = len(self.genes_covered.values())
        rep += "{0:d} gene(s) covered in {1:d} reads\n".format(ngenes, nreads)

        region_str = "\t{0:>9s} {1:7d}  {2:5.2f}%\n"
        rep += "Regions:\n"
        
        rep += region_str.format("5'-UTR", self.region_5prime,
                                 float(self.region_5prime)/float(self.reads_mapped)*100.0)
        rep += region_str.format("3'-UTR", self.region_3prime,
                                 float(self.region_3prime)/float(self.reads_mapped)*100.0)
        rep += region_str.format("coding", self.region_coding,
                                 float(self.region_coding)/float(self.reads_mapped)*100.0)
        rep += region_str.format("noncoding", self.region_noncoding,
                                 float(self.region_noncoding)/float(self.reads_mapped)*100.0)
        
        return rep
    pass

# For command-line options...
class Options:
    def __init__(self):
        self.gene_name = None
        self.match_limit = None
        self.show_multi_locus = False
    pass


class GeneNotFound(exceptions.Exception):
    def __init__(self):
        return
		
    def __str__(self):
        print "","An Gene Not Found exception occured!"

# main() takes an optional 'argv' argument, which allows us to call it
# from the interactive Python promp.

def main(argv = None):
    if argv is None:
        argv = sys.argv
    
    options = Options()
    stats = Stats()
    
    try:
        opts, args = getopt.getopt(argv, "hdml:g:", ["help", "output="])
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
        elif opt == '-g':  # limit the match to one particular gene
                           # takes a gene name as argument;
                           # -g NM_29391293
            options.gene_name = arg
        elif opt == '-l':  # limit the number of matching reads
                           # takes an integrer as argument
                           # stop processing after matching n reads to exons.
            options.match_limit = int(arg)
        elif opt == '-m':  # limit the number of matching reads
                           # takes an integrer as argument
                           # stop processing after matching n reads to exons.
            options.show_multi_locus = True
        else:
            usage()                         
            return 2

    if len(args) != 2:
        usage()                         
        return 2

    gene_file = args[0]
    reads_file = args[1]

    knowngenes = read_knowngenes(gene_file, options, stats)
    process_reads(reads_file, knowngenes, options, stats)

    #
    # print statistics.  Print to stderr so they do not interfere with
    # the normal output.
    #
    print >> sys.stderr, '\nStatistics:\n{0:s}'.format(stats)

    return 0




# In addition to reading in the exons for every known gene, this
# routine calculates two values that are used later.  The first is the
# position of the start of the exon w.r.t. the start of the gene.
# Exon #0 by definition starts at position 0.  
#
# The second is the position of the exon w.r.t. the start of the open
# reading from of the gene (aka CDS).  This value helps us calculate
# whether the ribosome falls on a codon boundary.
#

def read_knowngenes(known_genes, opt, stats):
    genename = opt.gene_name
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

# At one point I did not count non-coding genes.  Theoretically we
# should not have any reads from inside non-coding RNA because
# ribosomes don't read non-coding mRNA, but we seem to have hits
# there.
#
# If you really wanted to leave ncRNA out, you would have to filter
# the list of exons that are used to find overlaps in an earlier step
# to eliminate those that code for ncRNA.  That way you wouldn't see
# any overlap with ncRNA so the offending reads would be filtered out
# before this.
#
#
#             if gene.cds.start == gene.cds.end:
#                 # it is a non-coding RNA.
#                 # don't add to the gene list.
#                 continue
#            print "%s" % (g_row)
            
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
            gene.len = 0
            for i in xrange(int(g_row["exonCount"])):
                exon = Exon()
                exon.start = int(estart_list[i])
                exon.end = int(eend_list[i])
                exon.frame = int(efrm_list[i])
                exon_list.append(exon)
                gene.len += exon.len()

            # if reverse strand, reverse list of exons
            if  gene.strand == -1:
                exon_list.reverse()
            gene.exons = exon_list


            # calculate how far from the start of the gene each exon
            # begins.  The first exon in the correct direction will
            # start at 0.  We will use this later when calculating
            # where the ribosome is w.r.t. start of the transcript.
            #
            pos_gene = 0;
            for exon in exon_list:
                exon.pos_gene = pos_gene
                pos_gene += exon.len()
                    
            # calculate how far from the start of the coding
            # region each exon begins.
            #
            gene.annotate_cds()

            if _debug:
                print gene
                print ""

            # Some gene names appear in the genome more than once,
            # e.g. NM_001110250.  These entries have the same 
            # name but appear at different loci in the genome.  They
            # might have different introns but they should presumably
            # have the same sequence, or at least produce the same
            # amino acid sequence.
            #
            # These rare genes are the reason that each entry in this
            # gene dictionary (indexed by name) is a list.  Almost all
            # gene entries will be a list of exactly one entry, but a
            # few will have multiple entries in the list for a given
            # name.  Later when matching a read to a gene/exon name it
            # will be important to check which gene replicant the
            # read actually maps to.


            chrom = g_row["chrom"]
            if chrom not in chromosomes:
                chromosomes[chrom]= dict()
            if (gene.name not in chromosomes[chrom]):
                (chromosomes[chrom])[gene.name] = list()
            (chromosomes[chrom])[gene.name].append(gene)
            genecount += 1

            
    finally:
        if in_handle is not None:
            in_handle.close()

    return chromosomes
    # END read_knownegens()


# Given a name of a gene/exon pair, find the corresponding internal
# structure
def find_gene(knowngenes, row, opt, stats):

    ename = row['ename']
    # ename = NM_001011874_exon_2_0_chr1_3660633_r
    gname = ename[:ename.index('_exon')]
    # gname = NM_001011874

    genename = opt.gene_name
    if genename != None and genename != gname:
        raise GeneNotFound

    nexon = int( ename.split('_')[3])
    # enum = 2
    chrom = row['rchr']
    # chrom = chr1
    if chrom not in knowngenes:
        if chrom not in warnings:
            print >> sys.stderr, 'Warning:   Could not find known genes for chromosome "%s"' % (chrom)
            warnings[chrom] = True
            raise GeneNotFound

    # Some genes appear in the genome more than once,
    # e.g. NM_001110250.  We want to flag when that occurs.
    #
    # It is similar to when a read maps to multiple loci on
    # the genome, or maps to a single loci with multple genes.
    # It can be ambiguous which version of the gene the read
    # should be assigned to.  The default behavior is to skip
    # those reads but we give you the option to include them.
    #
    # We also collect and display statistics related to how
    # many reads mapped to genes like this that have multiple
    # loci.

    if (len((knowngenes[chrom])[gname]) > 1):
        stats.multi_locus_read += 1
        if (not opt.show_multi_locus):
            raise GeneNotFound

    rpos =  int(row['rpos'])
    rend = int(row['rend'])
    found = False
    for gene in (knowngenes[chrom])[gname]:
        if (gene.strand == -1):
            exon = gene.exons[len(gene.exons) - 1 - nexon]
        else:
            exon = gene.exons[nexon]
        if (exon.start <= rpos and exon.end >= rend):
            found = True
            break

    if (not found):
        print >> sys.stderr, 'Warning:   Could not find gene for read "%s""' % (row['rname'])
        raise GeneNotFound

    if (not gname in stats.genes_covered):
        stats.genes_covered[gname] = 1
    else:
        stats.genes_covered[gname] += 1
    
    return tuple([gene, exon])
    # END find_gene()




# Now read a file that is the output of intersectBed.
# It looks like this:
# chr1	3661187	3661209	ILLUMINA-D43DB5:7:10:17415:2809:0#0	255	-	chr1	3660632	3661579	NM_001011874_exon_2_0_chr1_3660633_r	0	-
#
# The first few fields are the read and the last few are the exon that
# this read maps to. 'junk' indicates junk field that we don't use.
# Symbolically, these might might be labelled:

# ---------- read ----------------- ------------ exon ---------------
# rchr rpos rend rname junk rstrand xchr xpos xend xname junk xstrand
#
# Should be able to get CDS of read.
# Scale that w/r/t gene length
# plot.

# process the reads
#
def process_reads(reads_file, knowngenes, opt, stats):
    match_limit = opt.match_limit

    hdr = "# {0:^35s} {1:^15s} {2:^8s} {3:^5s} {4:^5s} {5:^5s}"
    rec = "{0:37s} {1:15s} {2:8s} {3:5d} {4:5d} {5:6s} {6:5.2f}"

    try:
        warnings = dict()
        nmatched = 0

        in_handle = open(reads_file)

        print hdr.format("read name", "refseq", "gene", "dCDS", "dSOG", \
                         "regio", "%")
        
        reader = csv.DictReader(in_handle, delimiter='\t',
                                  fieldnames=['rchr', 'rpos', 'rend', 'rname', 'junk',
                                              'strand', 'echr', 'epos', 'eend', 'ename'])
        for row in reader:
            
            stats.reads_processed += 1

            try:
                gene, exon = find_gene(knowngenes, row, opt, stats)
            except GeneNotFound as (e):
                continue

            stats.reads_mapped += 1

            # Calculate distance from start of CDS.
            #
            dist_cdr_start = exon.cds.start
            if (gene.strand == 1):	# forward strand
                dist_cdr_start += int(row['rpos']) - exon.start
            else:			# reverse strand
                dist_cdr_start += exon.end - int(row['rend'])

            # Calculate distance from end of CDS.
            #
            dist_cds_end = exon.cds.end
            if (gene.strand == 1):	# forward strand
                dist_cds_end += int(row['rpos']) - exon.start
            else:			# reverse strand
                dist_cds_end += exon.end - int(row['rend'])

            if (dist_cds_end > 0):
                # length of 3'-UTR ...
                pos_pct = (float(dist_cds_end) / float(gene.len_utr_3())) * 100.0
                pos_which = "3'-UTR"
                stats.region_3prime += 1
            elif (dist_cdr_start < 0):
                # length of 5'-UTR ...
                len_5utr = gene.len_utr_5()
                # distance from start of gene ...
                dist_sog = len_5utr  - abs(dist_cdr_start)
                pos_pct = (float(dist_sog) / float(len_5utr)) * 100.0
                pos_which = "5'-UTR"
                stats.region_5prime += 1
            elif (gene.len_cds() > 0):
                pos_pct = float(dist_cdr_start) / float(gene.len_cds()) * 100.0
                pos_which = "CDS"
                stats.region_coding += 1
            else:
                pos_pct = float(dist_cdr_start) / float(gene.len) * 100.0
                pos_which = "NC"
                stats.region_noncoding += 1

            # calculate position w.r.t. beginning of the
            # transcript as a fraction of entire gene transcript length.
            #
            gpos = exon.pos_gene
            if (gene.strand == 1):	# forward strand
                gpos += int(row['rpos']) - exon.start
            else:			# reverse strand
                gpos += exon.end - int(row['rend'])
            
            pospct = (float(gpos)/float(gene.len)) * 100.0

            print rec.format(row['rname'], gene.name, gene.common_name, \
                             dist_cdr_start, gpos, pos_which, pos_pct)
                
            nmatched += 1
            if match_limit != None and nmatched >= match_limit:
                break

    # It's sort of unbelieveable, but python does not seem to have a
    # default handler for SIGPIPE.  We have to explicitly handle it
    # ourselves in case the output is being piped to a command and
    # that command finishes before we do.
    #

    except IOError as (e):
        if e.errno == errno.EPIPE:
            # EPIPE error
            pass
    finally:
        if in_handle is not None:
            in_handle.close()

    # END process_reads()


def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, "[-h] [-d] [-g gene_name] [-l #reads] [-m] <reference_tbl.csv> <alignments.bam> "

   
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
