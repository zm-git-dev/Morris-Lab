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
# ${SCRIPTS}/plot_rpos_vs_region.sh <tophat_out/aligned_position_stats.txt
#
#input format looks like this:
#    read name       refseq              gene     dCDS region    %    gene%
# ILLUMINA-D43DB5  NM_001164655    9530053A07Rik  3083  CDS    39.80  39.16
#
# read_name is the unique name assigned to thie read by the illumina sequencer
#
# refseq is the refseq name for the gene that this read maps to.
#
# gene and the common name for the gene that this read maps to.
#
# dCDS is how far (in nt) the read lies from the start of the CDS region of this gene.
#      Note that this value may be negative if in the 5'-UTR
#
# region is which region of the gene this read lies.  
#        Choices are 3-UTR, CDS, 5-UTR, NC.  NC is used for genes that are non-coding.
#
# % is how far (as a %) this read lies from the start of the designated region 
#   That is, how far in percentage terms from the start of the 3'-UTR, CDS, 5'-UTR region
#   In the case of non-coding genes (NC region), '%' is the same as 'gene%'.
#
# gene% is how far (as a %) this read lies from the start of the entire gene.
#

TMPFILE=`mktemp`

function on_exit()
{
    # guard against an empty TMPFILE variable wiping out everything

    if [ -n "$TMPFILE" ]; then 
	rm -f "$TMPFILE"*
    fi
}

trap on_exit EXIT


# Split the input file into several smaller files.
#
# For each type of region, make a temp file that contains all the
# reads that mapped to that type of region.  The result is four files,
#
#       $TMPFILE.3'-UTR, $TMPFILE.5'-UTR, $TMPFILE.CDS, $TMPFILE.NC
#
# each of which contains all the reads that mapped to one specific
# region across multiple genes.
#
# We also extract just the gene% and reseq name for each entry and place it in another 
# temp file, $TMPFILE.
#
# So from this one line, at least five temp files are made: four files
# categorized by region with read position as a percentage of that
# region, and a fifth file that contains data from all regions with
# read position as a percentage of gene length.
#

awk '!/^#/ {print $6,$2 >"'$TMPFILE'."$5; print $7,$2 >"'$TMPFILE'.gene" } ' -


# Calculate gene (or region) coverage at each position across all
# genes.  For each percentage position (0-100%) across all genes, show
# how many genes have at least one read at that position.  Also
# calculate coverage across percentage positions of each region of the
# gene (e.g. 3'-UTR, CDS, 5'-UTR).
#
# Perhaps one of the things you might learn from this is, say, few
# genes have reads mapped to 90% from the start of the gene.  Or maybe
# you would see that most of the reads are clustered at about 20% from
# the start of all genes.  Or this might show that there is even
# coverage of the gene as a whole but coverage of the 5'-UTR is
# clustered at the beginning.
#
# To help clarify, if there were 10 reads that mapped to 10 different
# positions on one gene, I would get 10 entries, one for each
# position, and each entry would have a value of 1 to represent the
# single gene that had coverage at that position.
#
# If there were 10 reads that mapped to 10 different positions on 10
# different genes then I would still get only 10 entries, one for each
# position, and each entry would still have a value of 1 to represent
# the single gene that had coverage at that position.  It doesn't
# matter if each position is represented by a different gene or the
# same gene; all that matters is that there is one gene that has
# coverage at each position.
#
# If there were 10 reads that mapped to only 2 positions on one gene,
# I would get 2 entries, one for each position, and each entry would
# have a value of one representing the single gene that had coverage
# at that position.
#
# If there were 10 reads that mapped to *the same position* on 5
# different genes, then I would get only one entry, but it would have
# a value of 5.
# 
# clear as mud, I know, but it is pretty much all done by this next command.

for i in $TMPFILE.*; 
do  
    sort -k 1n "$i" | uniq |  awk '{print $1}' | uniq -c >"$i.coverage" 
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
plot  "$TMPFILE.gene.coverage" using 2:1 title '' pt 13 ps 1

set origin 0.50, 0.50   # upper-right
set xlabel "Position in 5' region"
set title "5'-UTR position vs. Gene count"
plot  "$TMPFILE.5'-UTR.coverage" using 2:1 title '' pt 13 ps 1

set origin 0.50, 0.0   # lower-right
set xlabel "Position in 3' region"
set title "3'-UTR position vs. Gene count"
plot  "$TMPFILE.3'-UTR.coverage" using 2:1 title '' pt 13 ps 1

set origin 0.0, 0.0   # lower-left
set xlabel "Position in CDS region"
set title "CDS position vs. Gene count"
plot  "$TMPFILE.CDS.coverage" using 2:1 title '' pt 13 ps 1

unset multiplot
set term pop


EOF
