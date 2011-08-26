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
# 
# $ gnuplot ribopos_vs_genecnd.pl <out.txt


TMPFILE=`mktemp`

function on_exit()
{
    rm -f "$TMPFILE"*
}

trap on_exit EXIT

awk '!/^#/ {print >"'$TMPFILE'."$5; print $7,$2 >"'$TMPFILE'" } ' -

sort -k 1n $TMPFILE | uniq | awk '{print $1}' | uniq  -c >$TMPFILE.data

for i in $TMPFILE.*; 
do  
    awk '!/^#/ {print $6,$2}' <$i | sort -k 1n | uniq |  awk '{print $1}' | uniq -c >$i.data 
done

X11TERM="set terminal X11  size 1150,900"
PSTERM="set terminal postscript landscape color lw 0 \"Helvetica\" 12;  \
        set out 'plots.eps'"

PNGTERM="set terminal png;  \
        set out 'plots.png'"

gnuplot -persist <<EOF

set term push; 
$X11TERM

set multiplot
set size .5, .5

set autoscale                          # scale axes automatically
unset log                              # remove any log-scaling
unset label			       # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set ylabel "Number of genes"
set pointsize 1


set origin 0.0, 0.50   # upper-left
set xlabel "Position on gene"
set title "Whole gene position vs. Gene count"
plot  "$TMPFILE.data" using 2:1 title '' pt 13 ps 1

set origin 0.50, 0.50   # upper-right
set xlabel "Position in 5' region"
set title "5'-UTR position vs. Gene count"
plot  "$TMPFILE.5'-UTR.data" using 2:1 title '' pt 13 ps 1

set origin 0.50, 0.0   # lower-right
set xlabel "Position in 3' region"
set title "3'-UTR position vs. Gene count"
plot  "$TMPFILE.3'-UTR.data" using 2:1 title '' pt 13 ps 1

set origin 0.0, 0.0   # lower-left
set xlabel "Position in CDS region"
set title "CDS position vs. Gene count"
plot  "$TMPFILE.CDS.data" using 2:1 title '' pt 13 ps 1

unset multiplot
set term pop


EOF
