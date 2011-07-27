#!/bin/bash
# Workflow for aligning reads to the mm9 genome
#
# Usage:  flow1.sh <reads.fa>  
#

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

declare -A counts

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


log_msg "Trimming adapter and filtering bad reads..."

fastx_clipper -a ${ADAPTER} -l 20 -Q33 -v < ${rawreads} > ${basename}.multilen.fastq 2>clipstats.txt

# collect statistics
stats=($(awk '$2 ~ /^[0-9]+$/ {print $2 }' <clipstats.txt))
counts[raw]=${stats[0]}
counts[tooshort]=${stats[2]}
counts[adapter]=${stats[3]}
counts[uncalled]=${stats[4]}

log_msg "aligning to rRNA bowtie..."
bowtie --al contaminated.fq --un uncontaminated.fq /home/csw/Morris-Lab/ref-seq/contaminants/contaminants ${basename}.multilen.fastq -S  >${basename}.rrna.sam  

counts[rrna]=$(readc contaminated.fq)

# assert contaminated + uncontaminated == raw
reads=$(readc uncontaminated.fq)
log_msg "aligning ${reads} reads to  MM9 genome with tophat..."
tophat --segment-length 14 --segment-mismatches 0 -G ~/Morris-Lab/ref-seq/mm9/refseq_knowngenes.gtf  /mm9/mm9 uncontaminated.fq

# get how many failed to align
counts[unaligned]=$(( $(readc uncontaminated.fq) - $(samtools view tophat_out/accepted_hits.bam | awk '/\<NH:i:[0-9]+\>/ {print $1}' |sort|uniq|wc -l) ))

# how many reads mapped to multiple loci on the genome...
counts[nonuniq]=$(( $(samtools view tophat_out/accepted_hits.bam | awk '!/\<NH:i:1\>/ {print $1}' |sort|uniq|wc -l) ))

log_msg "filtering unique alignments..."
# Assemble only those alignments that were unique on the genome.
samtools view -h tophat_out/accepted_hits.bam | awk '/^@/ || /NH:i:1\>/'| samtools view -Sb -o tophat_out/accepted_hits_uniq.bam /dev/stdin 

log_msg "detecting exon converage..."
# intersect those alignmets with exons of known genes
intersectBed -abam tophat_out/accepted_hits_uniq.bam -b ${REFSEQ}/mm9/refseq_knowngenes_exons.bed -s -f 1.0 -wa -wb -bed > tophat_out/overlaps_exons.bed

# how many aligned reads fell outside of known exons
counts[introns]=$(( $(samtools view tophat_out/accepted_hits_uniq.bam | wc -l) - $(awk '{print $4}' <tophat_out/overlaps_exons.bed | sort|uniq|wc -l) ))

# count reads that overlap multiple features at a single locus
counts[mult_gene]=$(awk '{ print $4; }' <tophat_out/overlaps_exons.bed | sort | uniq -c | awk '$1 !~ /1/' | wc -l )

# Eliminate reads that map to more than a single gene
awk '{ print $_"\t"$4; }' <tophat_out/overlaps_exons.bed | sort -k 13 | uniq -f 12 -u | cut -f 1-12 >tophat_out/overlaps_exons_uniq.bed

log_msg "preparing coverage statistics..."
# prepare experimental statistics
python ${SCRIPTS}/rpos_dist.py ${REFSEQ}/mm9/refseq_knowngenes.txt tophat_out/overlaps_exons_uniq.bed >tophat_out/aligned_position_stats.txt

log_msg "graphing results..."
${SCRIPTS}/plot_rpos_vs_region.sh <tophat_out/aligned_position_stats.txt

printf "%20s\t%d\n" "Raw reads" ${counts[raw]}
printf "%20s\t%d\n" "uncalled base 'N'" ${counts[uncalled]}
printf "%20s\t%d\n" "adapter only" ${counts[adapter]}
printf "%20s\t%d\n" "too short (<20nt)" ${counts[tooshort]}
printf "%20s\t%d\n" "rRNA-spike" ${counts[rrna]}
printf "%20s\t%d\n" "unaligned" ${counts[unaligned]}
printf "%20s\t%d\n" "multiply aligned" ${counts[nonuniq]}
printf "%20s\t%d\n" "outside exons" ${counts[introns]}
printf "%20s\t%d\n" "multiple genes" ${counts[mult_gene]}

