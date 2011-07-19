# Plot the number of genes (y-axis) vs. read position the start of
# gene transcript (5' UTR) (x-axis)
#
# Usage:
#    ribopos_vs_genecount.sh <data.txt
#
# To prepare 'data.txt', start with a SAM file.  Typically output from
# tophat that has been filtered of reads that map to multiple
# locations (NM > 1).  Then run intersectBed to overlap those reads
# with exons to give overlaps_exons_uniq.bed Then run test.py to
# calculate the distance from start of transcript.
#
# $ python  scripts/test.py  ~/Morris-Lab/ref-seq/mm9/refseq_knowngenes.txt overlaps_exons_uniq.bed >out.txt
#
# w
# $ gnuplot ribopos_vs_genecnd.pl <out.txt


TMPFILE=`mktemp`
# TMPFILE=tmpfile

function on_exit()
{
    rm -f $TMPFILE
}

trap on_exit EXIT

awk '!/^#/ { print $7,$2}' - | sort -k 1n | uniq | awk '{print $1}' | uniq  -c >$TMPFILE 

gnuplot -persist <<EOF
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label			       # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "Position on gene"
set ylabel "Number of genes"
set title "Ribosome position vs. Gene count"


binwidth=5

bin(x,width)=width*floor(x/width)
#set xrange [0:100]


plot  "$TMPFILE" using 2:1 title '' 
#set bmargin 10 
#plot '$TMPFILE' using (bin(\$2,binwidth)):1 smooth freq title ''  with boxes

set term push
#set terminal postscript eps color lw 15 "Helvetica" 20
set terminal jpeg
set out 'a.jpg'
replot
set term pop


#plot '$TMPFILE'  using 2 title 'foo' 
#

EOF
