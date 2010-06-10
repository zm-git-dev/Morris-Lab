
import sys
import random
import getopt

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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
        sys.exit(2)                     
    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            usage()                     
            sys.exit()                  
        elif opt == '-d':                
            global _debug               
            _debug = 1                  
        elif opt in ("-o", "--output"):
            if arg=="":
                usage()                         
                sys.exit(2)                     
            output = arg
        elif opt == '-l':                
            mu = int(arg)

    if len(args) != 2:
        usage()                         
        sys.exit(2)                     

    if args[1].isdigit():
        count = int(args[1])
    else:
        usage()                         
        sys.exit(2)                     


    mkreads(args[0], output, count, mu, sigma)
    

def usage():
    print "usage: ", __name__, "[ -o <outfile> ] [ -l <read length> ] <reference genome>.fasta [count] usage message"
    
def mkreads(reference, outfile, count, mu, sigma):
    handle = open(reference, "rU")
    seqDict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()

    
    frags=[]
    
    sKeys = seqDict.keys()
    nSeq = len(sKeys)
    for i in range(count):
        # pick a sequence
        n = random.randint(0,nSeq-1)
        seq = seqDict[sKeys[n]]
        
        # pick a length of read to produce
        readlen = int(mu + np.random.randn()*sigma)

        # produce the read
        nBases = len(seq.seq)
        start = random.randint(0, nBases-readlen)
        read_frag = seq.seq[start:start+readlen]

        record = SeqRecord(read_frag,'fragment_%i' % (i+1),'','')
        frags.append(record)

    if outfile != "":
        output_handle = open(outfile, "w")
    else:
        output_handle = sys.stdout
    SeqIO.write(frags, output_handle, "fasta")
    if outfile != "":
        output_handle.close()


if __name__ == "__main__":
    main(sys.argv[1:])
