# compare_scatter.py
#
# Compare two datasets by making a scatterplot of both.
# This routine take a microsarray dataset from GEO and a rna-seq dataset
# and makea a scatter plot of the grene expression levels from each set.
#
# usage: python ../../scripts/compare_scatter.py platform.tbl refseq_reflink.tbl sample.tbl 
#
# platform.tbl - platform decsription file for microarray.
#		 e.g. http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL4045
#
# refseq_reflink.tbl - linkage file that allows us to map unigenes to refseq genes.
#		       We need this to map the unigene identifiers
#		       used by the microarray platform to something
#		       common like refseq genes in order to tell which
#		       genes we have in common.
#
# sample.tbl - the microarraye sample file.
#	       e.g. http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM476311
#
# rna-seq.tbl - the aligned short rna sequencing results from cufflinks.
#		e.g. /home/csw/Morris-Lab/data/mus-pro-11-30-2010/tophat_out/cufflinks_refseq/genes.expr
#


import csv
import sys
import getopt
import math
import string

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker

import numpy as np
from scipy import stats
import pylab
from array import array

# main() takes an optional 'argv' argument, which allows us to call it
# from the interactive Python promp.

def main(argv = None):
    if argv is None:
        argv = sys.argv

    xlabel = "X"
    ylabel = "Y"
    title = None
    
    try:
        opts, args = getopt.getopt(argv, "hx:y:d", ["help", "output="])
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
        elif opt == '-x':                
            xlabel = arg                  
        elif opt == '-y':                
            ylabel = arg                  
        elif opt == '-t':                
            title = arg                  
        else:
            usage()                         
            return 2


    # read in the list of x,y points for the plot from stdin
    #
    x = list()
    y = list()
    for line in sys.stdin.readlines():
        words = string.split(line) 
        x.append(float(words[0]))
        y.append(float(words[1]))

    global fig
    fig = plt.figure(1)
    ax = fig.add_subplot(111)


    # You might think that the way to get a log-log scatter plot
    # would be to use scatter() and set the axis scales to "log"
    # That doesn't work.  It ends up graphing only a portion of the data.
    # 	ax.set_yscale('log')
    #   ax.set_xscale('log')
    #   ax.set_xlim(1e-1, 10e5)
    #   ax.set_ylim(1e-1, 10e5)
    #   plt.scatter(x, y)
    #
    # Calling plt.loglog() works though.
    
    plt.loglog(x,y,marker='o',linestyle='none')

    # calculate regression statistics that we can display on the plot.
    #
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)


    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


    # Provide a default label if the user did not supply one.
    if (title == None):
        title = "Comparison of %s and %s" % (xlabel, ylabel)

    # Add a legend to display statistical measures.
    # Note that there is no guarantee that this will not overlap a datapoint.
    #
    ax.text(0.05, 0.9, 'R$\r{^{2}}$=%0.4f\nN=%d'%(r_value**2, len(x)),
            transform=ax.transAxes, va='top')
    
    ax.set_title(title)

    #fig.savefig("fig.%s.png" % (ylabel), format="png")
    
    plt.show()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
