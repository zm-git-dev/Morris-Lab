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

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker

import numpy as np
from array import array
from itertools import ifilter

import gnote
import genexref
import util



# main() takes an optional 'argv' argument, which allows us to call it
# from the interactive Python promp.

def main(argv = None):
    if argv is None:
        argv = sys.argv
    
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

    if len(args) != 4:
        usage()                         
        return 2



    array_list = prepare_microarray_list(args[0], args[1], args[2])
    seq_list = util.read_seq_list(args[3])

    global fig
    fig = plt.figure(1)
    scatterplot_adbundance(array_list, seq_list)
    plt.show()


    return 0



#     The purpose of this routine is to prepare a list of genes with
#     expression levels derived from microarray data.
#
#     The microarray sample files consist of lines like,
#
#     2 -0.4996132
#     3 0.07726434
#
#     This routine converts these to something like,
#
#        NM_133826	.707292
#        NM_177547	1.05501
#
#     The ID numbers has been replaced by unique gene names and the
#     values have been converted from log2 to linear.
#
#     We lookup each id in the platform file.  That will give use a
#     choice of gene annotations associated with that id.
#   

def prepare_microarray_list(platformFile, xrefFile, sampleFile):

    # The platform file allows us to map from microarray spot
    # identifier (a small integer) to a variety of gene annotations,
    # including various gene identifiers (unigene, refgene, etc) and
    # text annotation ("spermin binding protein")
    #
    platform_map = util.read_platform(platformFile, "UNIGENE")
    
    xref = genexref.GeneXref(xrefFile)

    abundance_dict = dict()
    sample_reader = csv.DictReader(open(sampleFile), delimiter='\t', fieldnames=['#id', '#abundance'])

    for row in sample_reader:
        try:
            entrez = platform_map[int(row['#id'])]
            gene = xref.unigene2refgene(entrez)
            value = math.pow(2, float(row['#abundance']))
            if gene in abundance_dict:
                abundance_dict[gene].append(value)
            else:
                abundance_dict[gene] = list([value])
                
            continue
        except ValueError:
            continue
        except KeyError:
            continue

    abundance_list = filter_by_stdev(abundance_dict)
    
    return abundance_list



def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <refseq.map> "





# 1 	0.1771207
# 2 	-0.4996132
# 3 	0.07726434
# 4 	-0.3684799
def build_sample_dict(sample_file, platform_map, xref):

    abundance_dict = dict()
    sample_reader = csv.DictReader(open(sample_file), delimiter='\t', fieldnames=['#id', '#abundance'])

    for row in sample_reader:
        try:
            entrez = platform_map[int(row['#id'])]
            gene = xref.unigene2refgene(entrez)
            value = math.pow(2, float(row['#abundance']))
            if gene in abundance_dict:
                abundance_dict[gene].append(value)
            else:
                abundance_dict[gene] = list([value])
                
            continue
        except ValueError:
            continue
        except KeyError:
            continue

    return abundance_dict

# Some of the entries in the microarray data are duplicates (at least
# according to their gene labels they are duplicates).  Sometimes
# these duplcates have highly variable expression levels and it is
# impossible to know which of the duplicates we should choose to
# represent a particular gene.

# What we do is based upon how many duplicate measurements there are for a particular gene.
# 	- If there are no duplicates, we are done.
# 	- If there are only two duplicates, the expression level is the average of the two duplciates.
# 	- If there are three or more duplicates, the expression levels
# 	  for identical genes are discarded entirely if the stddev is >
# 	  0.3.  If stddev < 0.3 then the mean is used.
#
# Input to this routine is a dictionary of tuples.  The key is a gene
# name.  Doesn't matter where this name comes from, as long as it is
# unique within the list of keys.  The value of each item in the
# dictionary is a variable length list of expression levels
# (linearized).
#
# This routine calculates the stddev and mean of these expression
# levels and filters the list down to one value.  If the stddev is <=
# 0.3 then we replace the list with a single value that is the mean of
# all the values for that item.  If the stddev is > 0.3 we toss the
# data entirely.
#
def filter_by_stdev(ad):
    out = list()
    for k in ad.keys():
        v = ad[k]
        if len(v) > 2:
            na = np.array(v)
            if na.std <= 0.3:
                out.append([k, na.mean()])
               #print "%s\t%d\t%f\t%f" % (k, len(v), na.std(), na.mean()) # equivalent to np.mean(x), np.var(x)

        elif (len(v) == 2):
            na = np.array(v)
            out.append([k, na.mean()])
            #print "%s\t%d\t%f\t%f" % (k, len(v), 1, na.mean()) # equivalent to np.mean(x), np.var(x)

        else:
            out.append([k, v[0]])
            #print "%s\t%d\t%f\t%f" % (k, len(v), 0, v[0]) # equivalent to np.mean(x), np.

    return out

def scatterplot_adbundance(arraydata_list, seqdata_list):
    """
    Make a histogram of normally distributed random numbers and plot the
    analytic PDF over it
    """
    ax = fig.add_subplot(111)

    a_dict = dict(arraydata_list)
    b_dict = dict(seqdata_list)

    # Given two sets of gene names, find which are common across both lists.
    # Those are the ones we want to put on the scatter plot.
    a = set(a_dict.keys())
    b = set(b_dict.keys())
    d = a.intersection(b)
    print "gene names in common = ", len(d)
    assert len(d)>0, 'No genes in common between array data and rna-seq data!'


    
    d_list = list(d)
    x = [a_dict[i] for i in d_list]
    y = [b_dict[i] for i in d_list]
    assert len(x)==len(y), 'length of arraydata and seqdata should be the same'
        
    plt.scatter(x, y)
    af =  gnote.AnnoteFinder(x,y, d_list)
    fig.canvas.mpl_connect('button_press_event', af)

    # Set the labels for the x-axis
    #ax.xaxis.set_major_locator(ticker.FixedLocator(x+.75))
    #ax.xaxis.set_major_formatter(ticker.FixedFormatter([item[0] for item in abundance_list]))

    # Set the angle of the X-axis labels.
    fig.autofmt_xdate(rotation=50)

    ax.set_xlabel('Nelson data')
    ax.set_ylabel('Morris data (FPKM)')
    ax.set_title(r'$\mathrm{Comparison\ of\ microarray\ data\ with\ RNA-Seq\ data}$')




if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
