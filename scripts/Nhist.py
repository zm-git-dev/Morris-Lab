#!python
# Count the number of times 'N' occurs at each position in a read.
# stdin should be just a list of reads, not a fasta or fastq file.
#
# usage:
# gunzip -c *A.gz | ~morrislab/scripts/reads | python ~morrislab/scripts/Nhist.py 082112_A


import sys
import re
import getopt


def main(argv = None):
    if argv is None:
        argv = sys.argv
    
    dataset_name = None
    
    try:
        opts, args = getopt.getopt(argv, "h", ["help", "output="])
    except getopt.GetoptError, msg:          
        usage(msg)                         
        return 2
    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            usage()                     
            sys.exit()                  
        elif opt == '--':  # end of options
            break
        else:
            usage()                         
            return 2

    dist = npos_distribution()
    if len(args) != 0:
        inject_npos_distribution(dist, args[0])

def usage(msg = None):
    if msg is not None:
        print msg
    print __name__, " <experiment_name> "

def npos_distribution():
    my_hist = {}

    for line in sys.stdin.readlines():
        starts = [match.start() for match in re.finditer("N", line)]
        for n in starts:
            if n in my_hist.keys():
                my_hist[n] += 1
            else:
                my_hist[n] = 1
    return my_hist

def inject_npos_distribution(dist, experiment):
    values = []
    for i in range(0,35):
        if i in dist.keys():
            values.append( dist[i])
        else:
            values.append( 0)
    print ', '.join(map(str, values))
   
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

 
