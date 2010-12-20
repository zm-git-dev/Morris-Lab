
import sys
import random
import getopt
import csv

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from array import array
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker
import numpy as np

def main(argv):
    mu = 28        # default average read length
    sigma = 5      # stdev of read length 
    count = 10     # default number of reads to produce
    output = ""    # output filename (defaults to stdout)
    try:
        opts, args = getopt.getopt(argv, "hl:o:d", ["help", "output="])
    except getopt.GetoptError:          
        usage()                         
        return 2                     
    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            usage()                     
            return 0
        elif opt == '-d':                
            global _debug               
            _debug = 1                  
        elif opt in ("-o", "--output"):
            if arg=="":
                usage()                         
                return 2
            output = arg
        elif opt == '-l':                
            mu = int(arg)

    if len(args) != 3:
        usage()                         
        return 2

    if args[2].isdigit():
        count = int(args[2])
    else:
        usage()                         
        return 2

    mkreads(args[0], args[1], output, count, mu, sigma)
    

def usage():
    print "usage: ", __name__, "[ -o <outfile> ] [ -l <read length> ] <reference genome>.fasta <reference_tbl.tsv> [count]"
    
def mkreads(reference, annotations, outfile, count, mu, sigma):
    x = array('l')

    an_dict = dict() # annotation dictionary
    frags=[]

    GeneReader = csv.DictReader(open(annotations), delimiter='\t')
    # this is bogus....  we need a way to randomly access this file.
    # If it is huge, we won't be able to read it into memory this way.
    # 
    for row in GeneReader:
        an_dict[row['#name']] = row
        

    # open the reference sequence file for random access
    fq_dict = SeqIO.index(reference, "fasta")
    nseq = len(fq_dict)


    for i in range(count):
        # choose a random sequence to read from.
        n = random.randint(0,nseq-1)
        seqid = fq_dict.keys()[n]
        seq = fq_dict[seqid].seq

        # choose a location within that sequence to read from. Start
        #          somewhere close to the 5' cds end.  Preferentially
        #          choose a position that is on the correct reading
        #          frame - specifically, a multiple of 3 nucleotides from
        #          the start of the coding region.


        # find the corresponding sequence in the annotation table.
        # That will tell us how long the sequence is and where the
        # coding region lies.
        seqinfo = an_dict[seqid]
        txStart = int(seqinfo['txStart'])
        cdsStart = int(seqinfo['cdsStart'])
        txStop = int(seqinfo['txEnd'])
        cdsStop = int(seqinfo['cdsEnd'])

        # pick a length of read to produce
        readlen = int(mu + np.random.randn()*sigma)
        
        # pick a start position close to the cdsstart
        startpos = int(np.random.randn()*5)

        # now pick the reading frame.  We want to try to get most of
        # the reads on a 3 codon boundary...
        
        if (startpos % 3) != 0:
            startpos = startpos + (3 - (startpos % 3))

        # gather some statistics on the distribution of start points.
        x.append(cdsStart - txStart)

        startpos += (cdsStart - txStart)
        if startpos < 0:
            startpos = 0
        
        endpos = startpos + readlen
        if endpos > txStop:
            endpos = txStop

        read_frag = seq[startpos:endpos]
        if len(read_frag) < readlen:
            xxx
        
        # record = SeqRecord(read_frag,'fragment_%i' % (i+1),'','')
        record = SeqRecord(read_frag,'%s_%i' % (fq_dict[seqid].name, i+1),'','')
        frags.append(record)


    print "avg = ", np. mean(x)
    print "stdev = ", np.std(x)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    n, bins, patches = plt.hist(x, 50, normed=0, facecolor='green', alpha=0.75)
    

    plt.show()
    
    if outfile != "":
        output_handle = open(outfile, "w")
    else:
        output_handle = sys.stdout
    SeqIO.write(frags, output_handle, "fasta")
    if outfile != "":
        output_handle.close()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
