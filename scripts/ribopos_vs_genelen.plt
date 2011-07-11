# To visualize the results from this routine,
#
# $ python  /home/csw/Morris-Lab/scripts/test.py  /home/csw/Morris-Lab/ref-seq/mm9/refseq_knowngenes.txt overlaps_exons_uniq.bed >out.txt
#
# $ sed '/^#/d' <out.txt | sort -k 6n | uniq -c -f 5 | awk '{print $1, $7}' >out2.txt
#
# $ gnuplot ribopos_vs_genelen.pl


set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label			       # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "Position on gene"
set ylabel "Number of reads"
set title "Read position vs. Read count"
set xr[40:60]
 set yr[0:1000]
plot  "out2.txt" using 2:1 title '' with lines
pause -1



