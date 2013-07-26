
# Convert a list of start positions and oligos into a GTF file
# suitable for use as an annottation file in IGV.  The output of this
# script can be used in the "Gene File" field when running "import
# genome" in IGV to make a new genome.

# The input to this script consists of a center position and an oligo, e.g.
# 12731 CACGTTCGTGTGGAACCTGGCGCTAAACCATTCGTAGACG
# 6975  GACACATTGATCATCGACACTTCGAACGCACTTGCGGCCC
# 8323  CGCTCGTGGGGGGCCCAAGTCCTTCTGATCGAGGCCCAGC
# 12268 AGCAGAATTCACCAAGCGTTGGATTGTTCACCCACTAATA
# 4891  TAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGG
# 5067  CGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGAC

{
    transcript++;
    start=$1-length($2)/2;
    end = start + length($2);
    printf "45S\t%d\t%d\toligo%d\n", start, end, transcript
}     

