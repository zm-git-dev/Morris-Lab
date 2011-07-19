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
    rm -f "$TMPFILE.*"
}

# trap on_exit EXIT

awk '{print >"'$TMPFILE'."$5 } ' -

for i in $TMPFILE.*; 
do  
    awk '!/^#/ {print $6,$2}' <$i | sort -k 1n | uniq |  awk '{print $1}' | uniq -c >$i.data 
done


gnuplot  <<EOF

# -persist
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label			       # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "Position in region"
set ylabel "Number of genes"


#set xrange [0:100]


set term x11 title "3'-UTR" position 0,100
set title "3'-UTR position vs. Gene count"
plot  "$TMPFILE.3'-UTR.data" using 2:1 title '' 

set term push
set terminal postscript landscape color lw 15 "Helvetica" 20
set size 1,1
set out '3primeUTR.eps'
replot
set term pop

set term x11 2 title "5'-UTR" position 1300,100
set title "5'-UTR position vs. Gene count"
plot  "$TMPFILE.5'-UTR.data" using 2:1 title '' 

set term push
set terminal postscript landscape color lw 15 "Helvetica" 20
set size 1,1
set out '5primeUTR.eps'
replot
set term pop

set term x11 1 title "CDS" position 650,100
set title "CDS position vs. Gene count"

binwidth=5
#set bmargin 10 
bin(x,width)=width*floor(x/width)
#plot "$TMPFILE.CDS.data" using (bin(\$2,binwidth)):1 smooth freq title ''  with boxes
set border -1 lw 0
plot  "$TMPFILE.CDS.data" using 2:1 title '' 

set term push
set terminal postscript landscape color lw 15 "Helvetica" 20
set size 1,1
set out 'CDS.eps'
replot
set term pop


EOF
