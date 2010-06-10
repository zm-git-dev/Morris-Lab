
# graph the frequency that reads align with genes.
#
# This routine will take as input a set of aligned reads and a data
# table made for the reference genome against which the reads were
# aligned.  The output is a histogram of frequencies that the reads
# aligned to each gene in the table.
#

import csv
import pysam
import sys
import getopt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
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

    if len(args) != 2:
        usage()                         
        return 2

    alignments, readcount = read_alignments(args[0], args[1])

    print "gene:\tfrequency"
    for a,b in readcount.iteritems():
        print "%s: %d" % (a, b)

    graph_frequency(alignments, readcount)
    return 0


def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <reference_tbl.csv> <alignments.bam> "




# Read the alignments and the reference genome.  Build a list of with
# one entry per aligned read that indictates which gene that read is
# aligned to.
#
def read_alignments(reference, alignFileBAM):

    readcount = dict()
    genes = list();
    alignments = array('l')

    
    GeneReader = csv.DictReader(open(reference), delimiter='\t')
    samfile = pysam.Samfile(alignFileBAM, "rb")

    for row in GeneReader:
        genes.append(row['#name'])
        readcount[row['#name']] = 0
        for alignedread in samfile.fetch( row['#name']):
            readcount[row['#name']] = readcount[row['#name']] + 1
            alignments.append(genes.index(row['#name']))

    return alignments, readcount


def graph_frequency(alignments, readcount):
    """
    Make a histogram of normally distributed random numbers and plot the
    analytic PDF over it
    """

    x = alignments


    fig = plt.figure()
    ax = fig.add_subplot(111)

    # the histogram of the data
    n, bins, patches = ax.hist(x, len(readcount), normed=1, facecolor='green', alpha=0.75)

    # hist uses np.histogram under the hood to create 'n' and 'bins'.
    # np.histogram returns the bin edges, so there will be 50 probability
    # density values in n, 51 bin edges in bins and 50 patches.  To get
    # everything lined up, we'll compute the bin centers
    #bincenters = 0.5*(bins[1:]+bins[:-1])
    # add a 'best fit' line for the normal PDF
    #y = mlab.normpdf( bincenters, mu, sigma)
    #l = ax.plot(bincenters, y, 'r--', linewidth=1)

    ax.set_xlabel('Genes')
    ax.set_ylabel('Frequency')
    #ax.set_title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
    ax.set_xlim(0, len(readcount))
    ax.set_ylim(0, .15)
    ax.grid(True)
    print "graphing ", len(alignments), " reads mapped onto ", len(readcount), " genes."

    plt.show()

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
