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
    nbins = 100
    regions = [ "5'-UTR", "CDS", "3'-UTR"  ]
    
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
        elif opt == '-r':                
            regions = arg.split()
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

    regionmaps = dict()
    for r in regions: regionmaps[r] = dict()
    in_handle = None
    try:
        in_handle = open(infile)

        binwidth = 100.0/nbins
        y = list()
        g_reader = csv.DictReader(in_handle, delimiter='\t')
        for row in g_reader:
            if row["regio"] in regions:
                bin = int(float(row['%']) / binwidth)

                regionmap = regionmaps[row["regio"]]
                if bin in regionmap:
                    regionmap[bin] += 1
                else:
                    regionmap[bin] = 1

    finally:
        if in_handle is not None:
            in_handle.close()


    fig = plt.figure(1)

    # Adjust space around the outside of the plots.
    fig.subplots_adjust(hspace=0.0001, wspace=0.0001,
                        bottom=0.07, 
                        left=0.05, right=0.96)
    
    # Provide a default title if the user did not supply one.
    if title is not None:
        fig.suptitle(title, fontsize=12)

    plotnum = 1
    for region in regions:
        regionmap = regionmaps[region]

        x = map(lambda k: (k*binwidth) + binwidth/2.0, regionmap.keys())
        y = map(lambda k: regionmap[k], regionmap.keys())

        if plotnum == 1:
            ax = fig.add_subplot(1,3,plotnum)
            ax1 = ax
        else:
            ax = fig.add_subplot(1,3,plotnum, sharey=ax1)

        plt.scatter(x,y,s=20,color='tomato')

        title = region
        ax.set_title(title)
        plt.xlim(0,100)

        if plotnum > 1:
            plt.setp(ax.get_yticklabels(), visible=False)
            plt.xticks((0,20,40,60,80,100),("","20","40","60","80","100"))
            
        xlabel = "Percent of %s region" % region
        ax.set_xlabel(xlabel)

        if plotnum is 1:
            ylabel = "# genes with mapped fragments"
            ax.set_ylabel(ylabel)
        plotnum += 1

    # show the command-line that invoked this plot.
    plt.figtext(.5,0.005, "  ".join(sys.argv), alpha=.3, ha='center')

    plt.show()

        
def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, "[-x <xlabel>] [-y <ylabel>] [-t <title>] [-m <marker>] [-n <#points>] <input_1.fpkm> <input_2.fpkm> "



if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
