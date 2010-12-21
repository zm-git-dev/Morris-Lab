# convert microarray data from unigene to refseq gene names and graph the expression level
#
# usage: python ../../scripts/readmap.py platform.tbl refseq_reflink.tbl sample.tbl  
# platform.tbl - platform tabel doenloaded from 


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

    if len(args) != 4:
        usage()                         
        return 2


    array_list = prepare_microarray_list(args[0], args[1], args[2])
    seq_list = read_seq_list(args[3])

    scatterplot_adbundance(array_list, seq_list)

    return 0


def prepare_microarray_list(platformFile, xrefFile, sampleFile):
    platform_map = read_map(platformFile)
    
    ref_map = read_unigene2refseq(xrefFile)

    ad = build_sample_dict(sampleFile, platform_map, ref_map)
    abundance_list = filter_by_stdev(ad)
    
    #     abundance_list = process_sample_file(args[2], platform_map, ref_map)
    #out = trim_abundance(abundance_list, -20);


    return abundance_list



def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <refseq.map> "




# Read the platform map
# 1	DV038825	243377
# 3	BM877548	76748
# 6	BC057190	104625
#
def read_map(reference): 
    platform = dict()
    genemap = csv.DictReader(open(reference), delimiter='\t',
                             fieldnames=[
                                 'ID', 'MetaCol', 'MetaRow', 'Column', 'Row', 'Biosequence Type',
                                 'Strand Type', 'Region', 'GB_ACC', 'SPOT_ID', 'UNIGENE',
                                 'Reporter Usage', 'Control Type', 'MPEDB Spot ID', 'Designation',
                                 'Related Gene Symbol', 'Description', 'GENE', 'Chromosome',
                                 'Name'])
    for row in genemap:
        try:
            platform[int(row['ID'])] = row['GENE']
        except ValueError:
            continue

    return platform


# Read the refseq_reflink.txt
#name	product	mrnaAcc	protAcc	geneName	prodName	locusLinkId	omimId
# MIR375		NR_035777		390601	0	100316383	0
# Os08g0544800	hypothetical protein	NM_001068944	NP_001062409	235443	97	4346214	0
# Os04g0476700	hypothetical protein	NM_001059617	NP_001053082	265619	97	4336157	0
# Os01g0170300	hypothetical protein	NM_001048674	NP_001042139	233119	97	4327633	0
#
def read_unigene2refseq(reflink_file):
    refmap = dict()
    refmap_reader = csv.DictReader(open(reflink_file), delimiter='\t')
    for row in refmap_reader:
        refmap[row['locusLinkId']] = row['mrnaAcc']

    return refmap


# 1 	0.1771207
# 2 	-0.4996132
# 3 	0.07726434
# 4 	-0.3684799
def process_sample_file(sample_file, platform_map, ref_map):

    abundance_list = list()
    sample_reader = csv.DictReader(open(sample_file), delimiter='\t', fieldnames=['#id', '#abundance'])

    for row in sample_reader:
        try:
            entrez = platform_map[int(row['#id'])]
            gene = ref_map[entrez]
            value = math.pow(2, float(row['#abundance']))
            abundance_list.append([gene, value])
            continue
        except ValueError:
            continue
        except KeyError:
            continue

    return abundance_list
                

# given a list like [('geneid', .0202) (geneid2, .1102)]
# filter the list.
#
# you can choose the top x percent or the bottom x percent,
# or choose just the top x or bottom x.
# if -1 < x < 1, x is a percentage.
# if x < -1 or x > 1 it is an absolute value.
#
# the new list is returned.
#
def trim_abundance(abundance_list, criteria):

    # make a copy of the input list because we are going to sort it
    # and I don't want to screw up the original.
    #
    genelist = abundance_list

    genelist.sort(key=lambda gene: gene[1])
    if (criteria > -1 and criteria < 1):
        # criteria is a percentage
        howmany = int(len(genelist) * criteria)
    else:
        # criteria is an absolute number of items
        howmany = criteria

    # if howmany < 0  we want to take the *last* howmany entries.
    # if howmany > 0 we want to take the *first* howmany entries.
    s = slice(howmany, None) if howmany < 0 else slice(0, howmany)

    # return a new sorted list containing just the slice indicated.
    return genelist.__getitem__(s)


def bargraph_adbundance(abundance_list):
    """
    Make a histogram of normally distributed random numbers and plot the
    analytic PDF over it
    """


    fig = plt.figure()
    ax = fig.add_subplot(111)

    x =  np.arange(len(abundance_list))
    l = [item[1] for item in abundance_list]
 
    plt.bar(x, l)

    # Set the labels for the x-axis
    ax.xaxis.set_major_locator(ticker.FixedLocator(x+.75))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter([item[0] for item in abundance_list]))

    # Set the angle of the X-axis labels.
    fig.autofmt_xdate(rotation=50)

    ax.set_xlabel('Genes')
    ax.set_ylabel('number of alignments')
    ax.set_title(r'$\mathrm{Count\ of\ read\ alignments\ per\ gene}$')

    plt.show()


# 1 	0.1771207
# 2 	-0.4996132
# 3 	0.07726434
# 4 	-0.3684799
def build_sample_dict(sample_file, platform_map, ref_map):

    abundance_dict = dict()
    sample_reader = csv.DictReader(open(sample_file), delimiter='\t', fieldnames=['#id', '#abundance'])

    for row in sample_reader:
        try:
            entrez = platform_map[int(row['#id'])]
            gene = ref_map[entrez]
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

#
# Read sequencing data
# 
# This should be the output of cufflinks run against the refseq knowgene database.
#
# gene_id	bundle_id	chr	left	right	FPKM	FPKM_conf_lo	FPKM_conf_hi	status
def read_seq_list(seq_file):

    seq_list = list()
    seq_reader = csv.DictReader(open(seq_file), delimiter='\t')

    for row in seq_reader:
        try:
            if (float(row['FPKM']) > 0):
                seq_list.append([row['gene_id'], float(row['FPKM'])])
            continue
        except ValueError:
            continue

    return seq_list


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
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    a = set([item[0] for item in arraydata_list])
    b = set([item[0] for item in seqdata_list])
    d = a.intersection(b)
    print "gene names in common = ", len(d)
    
    alist = filter(lambda item: item[0] in d,  arraydata_list)
    blist = filter(lambda item: item[0] in d,  seqdata_list)
        
    x = [item[1] for item in alist]     
    y = [item[1] for item in blist]     

    plt.scatter(x, y)

    # Set the labels for the x-axis
    #ax.xaxis.set_major_locator(ticker.FixedLocator(x+.75))
    #ax.xaxis.set_major_formatter(ticker.FixedFormatter([item[0] for item in abundance_list]))

    # Set the angle of the X-axis labels.
    fig.autofmt_xdate(rotation=50)

    ax.set_xlabel('Nelson data')
    ax.set_ylabel('Morris data (FPKM)')
    ax.set_title(r'$\mathrm{Comparison\ of\ microarray\ data\ with\ RNA-Seq\ data}$')

    plt.show()




if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
