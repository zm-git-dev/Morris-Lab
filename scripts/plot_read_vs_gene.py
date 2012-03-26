#!/usr/bin/python

import csv
import sys
import getopt
import string

import matplotlib.pyplot as plt

def main(argv = None):
    if argv is None:
        argv = sys.argv

    xlabel = None
    ylabel = None
    title = None
    marker=','
    nbins = None
    
    try:
        opts, args = getopt.getopt(argv, "hx:y:t:b:m:r:d", ["help", "output="])
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
        elif opt == '-b':                
            nbins = int(arg)
        elif opt == '-m':   # choose plot marker
            if arg.isdigit():
                marker = int(arg)
            else:
                marker = arg
        else:
            usage()                         
            return 2

    infile = args[0]

#input format looks like this:
#    read name       refseq              gene     dCDS region    %    gene%
# ILLUMINA-D43DB5  NM_001164655    9530053A07Rik  3083  CDS    39.80  39.16
#

# Calculate gene (or region) coverage at each position across all
# genes.  For each percentage position (0-100%) across all genes, show
# how many unique genes have at least one read at that position.  Also
# calculate coverage across percentage positions of each region of the
# gene (e.g. 3'-UTR, CDS, 5'-UTR).

    regionmap = dict()
    in_handle = None
    try:
        in_handle = open(infile)

        g_reader = csv.DictReader(in_handle, delimiter='\t')
        for row in g_reader:
            if nbins is not None:
                bin = float(row['gene%']) / (100.0/nbins)
                bin = (int(bin) * (100.0/nbins)) + (100.0/nbins)/2.0
            else:
                bin = float(row['gene%'])

            if bin in regionmap:
                regionmap[bin] += 1
            else:
                regionmap[bin] = 1

    finally:
        if in_handle is not None:
            in_handle.close()


    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    # Adjust space around the outside of the plots.
    fig.subplots_adjust(hspace=0.0001, wspace=0.0001,
                        bottom=0.07, 
                        left=0.05, right=0.96)
    

    # Provide a default title if the user did not supply one.
    if title is not None:
        fig.suptitle(title, fontsize=12)

    x = list()
    y = list()
    for k in regionmap.keys():
        x.append(k)
        y.append(regionmap[k])

    plt.scatter(x,y,s=20,color='tomato')

    plt.xlim(0,100)

    ax.set_xlabel("Percent of gene")

    ax.set_ylabel("# reads")

    plt.show()

        
def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, "[-x <xlabel>] [-y <ylabel>] [-t <title>] [-m <marker>] [-n <#points>] <input_1.fpkm> <input_2.fpkm> "



if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
