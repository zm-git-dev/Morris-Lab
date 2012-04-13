#!/usr/bin/python
# fpkm_scatter.py
# 
#
# Compare two datasets by making a scatterplot of both.
# This is an incredibly unsophisticated script that just creates a
# scatter plot from two sets of data.  No complicated processing or
# filtering is happening here.
#
# The input files must be text files, each having a number of columns
# seperated by a single tab, and the first row must be a header with
# the names of the columns.  Two of the columns must be labelled
# 'tracking_id' and 'FPKM'.  Other columns will be ignored.
#
# The following would be a valid format for the input files,
#   tracking _id  col2    col3    FPKM
#   gene1         ajunk   bjunk   10.292
#   gene2         junkk   jjunk   0.2e-12
#
# The columns must be labelled, seperated by tabs, and may be in any
# order, and may include junk columns.  
#
# usage: python ../../scripts/cuffscatter.py input_A input_B
#
# The input_A and input_B may be the output of cufflinks
# ('genes.fpk.tracking') or they may be the output of more stringent
# processing by 'rpos_dist.py' followed by 'calculate_fpkm.pk'
#
# Options:
#
# -x "label"	Use 'label' for the X axis label; default is no label
# -y "label"	Use 'label' for the Y axis label; default is no label
# -t "title"	Use 'title' for the title of the plot; default is none,
#		or if -x and -y are supplied, then a title will be constructed
#		from the X- and Y- labels.
# -n <num>	Graph only the top <num> most abundant genes.	
#
# -m <marker>	Mark each point with <marker>; default is a single pixel mark.
# 		<marker> can be any of the matplotlib line markers documented at
# 		http://matplotlib.sourceforge.net/api/artist_api.html#matplotlib.lines.Line2D.set_marker
#		Use -m 'o' for a larger round mark that is easier to see.
#


import csv
import sys
import getopt
import math
import copy
import string

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker

import numpy as np
from scipy import stats
import pylab
from array import array

class AnnoteFinder:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.
    
  Register this function like this:
    
  scatter(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, xdata, ydata, annotes, axis=None):
    self.data = zip(xdata, ydata, annotes)
    if axis is None:
      self.axis = pylab.gca()
    else:
      self.axis= axis

  def distance(self, x1, x2, y1, y2):
    """
    return the distance between two points
    """
    return math.hypot(x1 - x2, y1 - y2)

  def findAnnote(self, clickX, clickY):
      annote=None
      annotes = []
      for x,y,a in self.data:
             annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
      if annotes:
          annotes.sort()
          distance, x, y, annote = annotes[0]
          return annote
          # return '(%3.2f, %3.2f) (%3.2f, %3.2f) %s'%(x, y, clickX, clickY, annote)



#####################################################################
# Class Gene - Class used to hold data for one gene.  One of these is
# generated for each line of the input files.
#
class Gene:
    def __init__(self):
        self.name = "name"
        self.common_name = "common_name"
        self.len = 0
        self.nfrags = 0
        self.coding = 0
        
    def __str__(self):
        rep = '{0} len: {2} \n'\
              .format(self.name, 
                      self.len )

        return rep

    # End Gene


#####################################################################
# main() takes an optional 'argv' argument, which allows us to call it
# from the interactive Python promp.
#
def main(argv = None):
    if argv is None:
        argv = sys.argv

    xlabel = None
    ylabel = None
    title = None
    nsamples = None
    include_missing = False
    marker=','
    
    try:
        opts, args = getopt.getopt(argv, "hx:y:t:n:m:dz", ["help", "output="])
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
        elif opt == '-n':                
            nsamples = int(arg)
        elif opt == '-z':                
            include_missing = True
        elif opt == '-m':
            if arg.isdigit():
                marker = int(arg)
            else:
                marker = arg
        else:
            usage()                         
            return 2


    cuffout1 = args[0]
    cuffout2 = args[1]


    # read in the list of x,y points for the plot from stdin
    #
    genes1 = dict()
    in_handle = None
    try:
        in_handle = open(cuffout1)
        
        g_reader = csv.DictReader(in_handle, delimiter='\t')
        for row in g_reader:
            gene = Gene()
            gene.formal_name = row["tracking_id"]
            gene.common_name = row["gene"]
            gene.fpkm = float(row["FPKM"])
            gene.coding = int(row["coding"])
            genes1[gene.formal_name] = gene

    finally:
        if in_handle is not None:
            in_handle.close()


    # read in the list of x,y points for the plot from stdin
    #
    genes2 = dict()
    in_handle = None
    try:
        in_handle = open(cuffout2)
        
        g_reader = csv.DictReader(in_handle, delimiter='\t')
        for row in g_reader:
            gene = Gene()
            gene.formal_name = row["tracking_id"]
            gene.common_name = row["gene"]
            gene.fpkm = float(row["FPKM"])
            gene.coding = int(row["coding"])
            genes2[gene.formal_name] = gene

    finally:
        if in_handle is not None:
            in_handle.close()

    # If the user wants, we can include the values that are missing in one dataset or another.
    #
    if include_missing:
        geneset1 = set(genes1.keys())
        geneset2 = set(genes2.keys())
        for k in geneset1 - geneset2:
            gene = copy.copy(genes1[k])
            gene.fpkm = 0
            genes2[k] = gene

        for k in geneset2 - geneset1:
            gene = copy.copy(genes2[k])
            gene.fpkm = 0
            genes1[k] = gene
    else:
        for k in set(genes1.keys()).difference(genes2.keys()):
            del genes1[k]

    keys = sorted(genes1, key=lambda key: genes1[key].fpkm+genes2[key].fpkm, reverse=True)

    if nsamples is not None:
        nsamples = min(nsamples, len(keys))
    else:
        nsamples = len(keys)
        
    # print keys[:nsamples]
    
    x = list()
    y = list()
    annotes = list()
    colors = list()
    for k in keys[:nsamples]:
        gene = genes1[k]
        x.append(gene.fpkm)
        y.append(genes2[k].fpkm)
        annotes.append(gene.common_name+'  '+gene.formal_name)
        if (gene.coding):
          colors.append('blue')
        else:
          colors.append('tomato')
        
    global fig
    fig = plt.figure(1)
    ax = fig.add_subplot(111)


    # You might think that the way to get a log-log scatter plot
    # would be to use scatter() and set the axis scales to "log"
    # That doesn't work.  It ends up graphing only a portion of the data.
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.scatter(x, y, c=colors)
    ax.autoscale_view(False, True, True)

    #
    # Calling plt.loglog() works though.
    
    #plt.loglog(x,y,marker=marker,linestyle='none', color='tomato')

    # calculate regression statistics that we can display on the plot.
    #
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
        
    if ylabel is not None:
        ax.set_ylabel(ylabel)


    # Add a legend to display statistical measures.
    # Note that there is no guarantee that this will not overlap a datapoint.
    #
    ax.text(0.05, 0.9, 'R$\r{^{2}}$=%0.4f\nN=%d'%(r_value**2, len(x)),
            transform=ax.transAxes, va='top')
    
    # Provide a default title if the user did not supply one.
    if title is None and xlabel is not None and ylabel is not None:
        title = "Comparison of %s and %s" % (xlabel, ylabel)
    if title is not None:
        ax.set_title(title)

    # Print the command-line that invoked this plot.
    plt.figtext(.5,0.005, " ".join(sys.argv), alpha=.3, ha='center')

    af =  AnnoteFinder(x,y, annotes)
    ax.format_coord = lambda x,y : af.findAnnote(x, y)

    #fig.savefig("fig.%s.png" % (ylabel), format="png")
    
    plt.show()

    return 0

def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, "[-x <xlabel>] [-y <ylabel>] [-t <title>] [-m <marker>] [-n <#points>] <input_1.fpkm> <input_2.fpkm> "


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
