
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
import matplotlib.ticker as ticker
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

    startpos = read_startpos(args[0], args[1])

    #print "start pos:\tfrequency"
    #for a,b in startpos.iteritems():
        #print "%s: %d" % (a, b)

    graph_startpos(startpos)
    return 0


def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <reference_tbl.csv> <alignments.bam> "




# Read the the start position of each alignment w.r.t. the coding
# reagion of the gene.  return a dictionary whose keys are relative
# start positions ('-89', '288', '-22', etc.) and whose values are the
# number of reads that we saw start at that position.
#
def read_startpos(reference, alignFileBAM):

    startpos = dict()
    
    GeneReader = csv.DictReader(open(reference), delimiter='\t')
    samfile = pysam.Samfile(alignFileBAM, "rb")

    for row in GeneReader:
        for alignedread in samfile.fetch( row['#name']):
            pos = alignedread.pos
            txStart = int(row['txStart'])
            cdsStart = int(row['cdsStart'])
            
            dStart = pos # (txStart + pos) - cdsStart
            
            #if dStart < -1000:
            #    print "%s\t%d\t%d\ttx=%d\tcds=%d\tUTR=%d" % \
            #          (row['#name'], dStart, pos, txStart, cdsStart, 
            #          (cdsStart - txStart))
            
            key = str(dStart)
            if key in startpos:
                startpos[key] += 1
            else:
                startpos[key] = 1

    return startpos


def graph_startpos(startpos):
    """
    Make a histogram of normally distributed random numbers and plot the
    analytic PDF over it
    """

    # find the minimum and maximum start positions that were seen.  We
    # will use this information for setting the x-axis range.
    e = [ int(x) for x in startpos.keys() ]
    minPos = min(e)
    maxPos = max(e)
    print minPos, maxPos

    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    x = np.arange(minPos, maxPos)

    y =  map(lambda k: startpos.get(str(k), 0), x)
    
    plt.plot(x, y)
    plt.xlim(minPos, maxPos)
    

    ax.set_xlabel('CDS Position [nt from start]')
    ax.set_ylabel('Read 5\' ends')
    ax.set_title(r'$\mathrm{Count\ of\ read\ alignments\ per\ gene}$')

    plt.show()

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
