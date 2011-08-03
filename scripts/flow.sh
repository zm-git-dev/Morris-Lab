#!/bin/bash
# Workflow for aligning reads to the mm9 genome
#
# Usage:  flow1.sh <reads.fa>  
#

usage()
{
cat << EOF
usage: $0 options

This script run the test1 or test2 over a machine.

OPTIONS:
   -h      Show this message
   -t      Test type, can be ‘test1′ or ‘test2′
   -r      Server address
   -p      Server root password
   -v      Verbose
EOF
}

VERBOSE=
STATONLY=
while getopts "hsv" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         s)
             STATSONLY=1
	     echo "stats only"
             ;;
         v)
             VERBOSE=1
             ;;
         ?)
             usage
             exit
             ;;
     esac
done


logfile=/tmp/mylog

echo >$logfile
trap "rm -f $logfile" EXIT


# trap control-C or else it will only interrupt the currently running
# command and the script will continue to run.
trap 'exit 0' INT

# Output message to log file.
function log_msg()
{
    echo "$*" >>$logfile
}

# Count reads in a FASTQ file
function readc() { echo $(( $(cat $* | wc -l) / 4 )); }

exec 2>errors.txt

declare -A reads
declare -A alignments

ADAPTER="ATCTCGTATGCCGTCTTCTGCTTG"
MORRIS="/home/csw/Morris-Lab"
REFSEQ="${MORRIS}/ref-seq"
SCRIPTS="${MORRIS}/scripts"

rawreads="$1"

# everything after last '/'
basename=${rawreads##*/}
# everything before last '.'
basename=${basename%.*}

# If output is to a terminal, Start spinner
if [ -t "1" ]; then
    ${SCRIPTS}/spinner.sh &
fi


if [[ "$STATSONLY" -ne 1 ]]
then
    log_msg "Trimming adapter and filtering bad reads..."

    fastx_clipper -a ${ADAPTER} -l 20 -Q33 -v < ${rawreads} > ${basename}.multilen.fastq 2>clipstats.txt
    echo "fastx_clipper returned $?" >&2

    log_msg "aligning to rRNA bowtie..."
    bowtie --al contaminated.fq --un uncontaminated.fq /home/csw/Morris-Lab/ref-seq/contaminants/contaminants ${basename}.multilen.fastq -S  >${basename}.rrna.sam  
    echo "bowtie returned $?" >&2

    goodreads=$(readc uncontaminated.fq)
    log_msg "aligning ${goodreads} reads to  MM9 genome with tophat..."
    tophat --segment-length 10 --segment-mismatches 0 -G ~/Morris-Lab/ref-seq/mm9/refseq_knowngenes.gtf  /mm9/mm9 uncontaminated.fq
    echo "tophat returned $?" >&2

    log_msg "filtering out imperfect alignments..."
    samtools view -h tophat_out/accepted_hits.bam | awk '/^@/ || /NM:i:0\>/'| samtools view -Sb -o tophat_out/accepted_hits_perfect.bam /dev/stdin 
    
    log_msg "detecting exon converage..."
    # intersect those alignmets with exons of known genes
    intersectBed -abam tophat_out/accepted_hits_perfect.bam -b ${REFSEQ}/mm9/refseq_knowngenes_exons.bed -s -f 1.0 -wa -wb -bed > tophat_out/overlaps_exons.bed
    echo "intersectBed returned $?" >&2

    # Eliminate reads that map to more than a single gene
    awk '{ print $_"\t"$4; }' <tophat_out/overlaps_exons.bed | sort -k 13 | uniq -f 12 -u | cut -f 1-12 >tophat_out/overlaps_exons_uniq.bed
fi


# collect statistics
stats=($(awk '$2 ~ /^[0-9]+$/ {print $2 }' <clipstats.txt))
reads[raw]=${stats[0]}
reads[tooshort]=$(( ${reads[raw]} - ${stats[2]} ))
reads[adapter]=$(( ${reads[tooshort]} - ${stats[3]} ))
reads[uncalled]=$(( ${reads[adapter]} - ${stats[4]} ))


log_msg "Counting rRNA alignments..."
reads[rrna]=$(readc uncontaminated.fq)

# Count how many failed to align
log_msg "Counting failed alignments..."
reads[unaligned]=$(( $(samtools view tophat_out/accepted_hits.bam | awk '{print $1}' |sort|uniq|wc -l) ))
alignments[unaligned]=$(samtools view tophat_out/accepted_hits.bam | wc -l )


# Count how many imperfect alignments
log_msg "Counting imperfect alignments..."
reads[imperfect]=$(($(samtools view tophat_out/accepted_hits_perfect.bam | awk '{print $1}' | sort|uniq|wc -l) ))

alignments[imperfect]=$(samtools view tophat_out/accepted_hits_perfect.bam | wc -l)

# how many aligned reads fell outside of known exons
log_msg "Counting reads outside of exons..."
reads[introns]=$(awk '{print $4}' <tophat_out/overlaps_exons.bed | sort|uniq|wc -l)
alignments[introns]=$(wc -l <tophat_out/overlaps_exons.bed)

# count reads that overlap multiple features at a single locus
log_msg "Counting overlaps with multiple features..."
reads[mult_gene]=$(awk '{ print $4; }' <tophat_out/overlaps_exons_uniq.bed | sort | uniq -c | wc -l )

alignments[mult_gene]=$(wc -l <tophat_out/overlaps_exons_uniq.bed)
    
log_msg "preparing coverage statistics..."
# prepare experimental statistics
python ${SCRIPTS}/rpos_dist.py ${REFSEQ}/mm9/refseq_knowngenes.txt tophat_out/overlaps_exons_uniq.bed  >tophat_out/aligned_position_stats.txt 2>rpos_dist_stats.txt
echo "python rpos_dist.py returned $?" >&2

log_msg "graphing results..."
${SCRIPTS}/plot_rpos_vs_region.sh <tophat_out/aligned_position_stats.txt
log_msg " "
sleep 3
echo 

printf "%20s\t%d\t%d\n" "Raw reads" ${reads[raw]} ${alignments[raw]}
printf "%20s\t%d\t%d\n" "too short (<20nt)" ${reads[tooshort]} ${alignments[tooshort]}
printf "%20s\t%d\t%d\n" "adapter only" ${reads[adapter]} ${alignments[adapter]}
printf "%20s\t%d\t%d\n" "uncalled base 'N'" ${reads[uncalled]} ${alignments[uncalled]}
printf "%20s\t%d\t%d\n" "rRNA-spike" ${reads[rrna]} ${alignments[rrna]}
printf "%20s\t%d\t%d\n" "unaligned" ${reads[unaligned]} ${alignments[unaligned]}
printf "%20s\t%d\t%d\n" "imperfect alignments" ${reads[imperfect]} ${alignments[imperfect]}
printf "%20s\t%d\t%d\n" "outside exons" ${reads[introns]} ${alignments[introns]}
printf "%20s\t%d\t%d\n" "multiple genes" ${reads[mult_gene]} ${alignments[mult_gene]}

echo
cat rpos_dist_stats.txt