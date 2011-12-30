#!/bin/bash
# Workflow for aligning reads to the mm9 genome
#
# Usage:  flow1.sh <reads.fa>  
#

usage()
{
cat << EOF

Align reads to reference genome and calculate position of read 5' ends
w.r.t. beginning of CDS region.  The input data is extensively
filtered to eliminate low quality reads and reads that map to
ribosomal RNA.

usage: `basename $0` [options] <data.fastq>

e.g.,

    `basename $0` liver.fastq 
    Fully process set of short reads from liver.


    `basename $0` -s -g yeast SRR014368.fastq 
    Only calculate statistics based upon a previous complete run.  Use
    s_cerevisiae reference genome.

OPTIONS:
   -h      Show this message
   -s      calculate statistics only
   -g      reference genome (mm9, hg18, or hg19)
   -m      number of mismatches to allow (default is zero)
   -b      use bowtie instead of tophat (not tested)
   -a      adapter (default is ATCTCGTATGCCGTCTTCTGCTTG)
EOF
}

command_line="$0 $@"

USE_BOWTIE=
STATONLY=
LENGTH=20
MISMATCHES=$(( 0 ))
REF_GENOME=mm9
ADAPTER="ATCTCGTATGCCGTCTTCTGCTTG"
while getopts "hsvbm:g:a:l:" OPTION
do
     case "${OPTION}" in
         h)
             usage
             exit 1
             ;;
         s)
             STATSONLY=1
	     echo "stats only"
             ;;
         g)
 	     REF_GENOME="${OPTARG}"
	     ;;
         m)
 	     MISMATCHES=$(( ${OPTARG} ))
	     ;;
         b)
 	     USE_BOWTIE=1
	     ;;
         a)
	     ADAPTER="${OPTARG}"
	     ;;
         l)
	     LENGTH="${OPTARG}"
	     ;;
         ?)
             usage
             exit
             ;;
     esac
done

case "${REF_GENOME}" in 
    "mm9")
 	INDEX_BASE="mm9/mm9"
	RRNA_BASE="mm9/rrna"
	REFSEQ_BASE="mm9"
	;;
    "hg18")
 	INDEX_BASE="hg18/hg18"
	RRNA_BASE="hg18/rrna"
	REFSEQ_BASE="hg18"
	;;
    "hg19")
 	INDEX_BASE="hg19/hg19"
	RRNA_BASE="hg19/rrna"
	REFSEQ_BASE="hg19"
	;;
    "yeast"|"s_cerevisiae")
 	INDEX_BASE="s_cerevisiae/s_cerevisiae"
	RRNA_BASE="s_cerevisiae/rrna"
	REFSEQ_BASE="s_cerevisiae"
	;;
    *)
	echo REFGENOME
	usage
	exit
	;;
esac

# getopts will not change the positional parameter set — if you want
# to shift it, you have to do it manually after processing:
shift $((OPTIND-1))


rawreads="$1"
echo "Raw reads is ${rawreads}"

#
# Use the name of the file that contains the raw reads to generate the
# base name of all our output files.
# 
basename=`basename ${rawreads}`
basename=${basename%.*}
echo "basename =:$basename:"

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

function filter_alignments()
{
    # "filtering out imperfect alignments..."
    samtools view -h /dev/stdin   \
	| awk '# filter out alignments with more than ${MISMATCHES} mismatched bases.
               # used to produce accepted_hits_perfect.bam
               /^@/ {print}
               !/^@/ { for (i = 1; i <= NF; i++) 
	                 if ($i ~ /\<NM:i:[0-9]+\>/ && match($i, /[0-9]+/) \
                             && substr($i, RSTART, RLENGTH) >= '"${MISMATCHES}"') {
 	                      print;
                              count++;
                              break;
                          }
                      }
		END {
		    print "reads with fewer than '"${MISMATCHES}"' mismatches:\t", count > "/dev/stderr" # gawk only
		    }' \
    	| awk '# filter out alignments with junctions...
               #  used to produce accepted_hits_nojunc.bam
               /^@/ { print }
               $6 ~ /^[0-9]+M$/ {
                      print;
                      count++;
                      }
               END {
                   print "reads without junctions:\t",count > "/dev/stderr"   # gawk only
                   }' \
	| samtools view -Sb -o /dev/stdout /dev/stdin \
	| intersectBed -abam /dev/stdin -b ${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes_exons.bed -s -f 1.0 -wa -wb -bed \
        | awk '# Eliminate reads that map to more than a single gene
               # used to produce overlaps_exons_uniq.bed
               { print $_"\t"$4; }' | sort -k 13 | uniq -f 12 -u | cut -f 1-12 

    # echo "intersectBed returned $?" >&2


    mismatches=$*
}

# Count reads in a FASTQ file
function readc() { echo $(( $(cat $* | wc -l) / 4 )); }

exec 2>errors.txt

declare -A reads
declare -A alignments

MORRIS="/home/csw/Morris-Lab"
REFDIR="/usr/local/share/genome"
SCRIPTS="${MORRIS}/scripts"

# If output is to a terminal, Start spinner
if [ -t "1" ]; then
    ${SCRIPTS}/spinner.sh &
fi

for i in ${rawreads} ${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes.gtf ${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes_exons.bed
do
    if [ ! -r "$i" ]
    then
	echo Cannot find input file "$i".
	echo This file is necessary to process this sample.
	exit 1;
    fi
done
       

if [[ "$STATSONLY" -ne 1 ]]
then
    echo "${command_line}" >flow.cmd

    log_msg "Trimming adapter, clipping to ${LENGTH}nt, and filtering bad reads..."

    fastx_clipper -a "${ADAPTER}" -l ${LENGTH} -Q33 -v < ${rawreads} > ${basename}.multilen.fastq 2>clipstats.txt
    echo "fastx_clipper returned $?" >&2

    log_msg "aligning to rRNA bowtie..."
    bowtie --al rrna.fq --un non-rrna.fq "${RRNA_BASE}" ${basename}.multilen.fastq -S  >${basename}.rrna.sam  
    echo "bowtie returned $?" >&2

    goodreads=$(readc non-rrna.fq)
    if [ -z "${USE_BOWTIE}" ]
    then
	log_msg "aligning ${goodreads} reads to ${REF_GENOME} genome with tophat..."
	tophat --segment-length 10 --segment-mismatches 0 -G ${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes.gtf  "/${INDEX_BASE}"  non-rrna.fq
	echo "tophat returned $?" >&2
    else
	echo "using bowtie for main aliagnment" >&2
	log_msg "aligning ${goodreads} reads to ${REF_GENOME} genome with bowtie..."
	bowtie --best --sam "/${INDEX_BASE}"  non-rrna.fq | awk '$3!="*"' >accepted_hits.sam
	echo "bowtie returned $?" >&2
	mkdir tophat_out
	samtools view -Sb -o tophat_out/accepted_hits.bam accepted_hits.sam
    fi

    log_msg "filtering out imperfect alignments..."
    samtools view -h tophat_out/accepted_hits.bam        \
         | awk '/^@/ {print}
                !/^@/ { for (i = 1; i <= NF; i++) 
	                  if ($i ~ /\<NM:i:[0-9]+\>/ && match($i, /[0-9]+/) && substr($i, RSTART, RLENGTH) >= '"${MISMATCHES}"') {
 	                      print ;
                              break;
                          }
                      }' \
	 | samtools view -Sb -o tophat_out/accepted_hits_perfect.bam /dev/stdin 
    
    log_msg "filtering out alignments with junctions..."
    samtools view -h tophat_out/accepted_hits_perfect.bam | awk '/^@/ || $6 ~ /^[0-9]+M$/'| samtools view -Sb -o tophat_out/accepted_hits_nojunc.bam /dev/stdin 
    
    log_msg "detecting exon converage..."
    # intersect those alignmets with exons of known genes
    intersectBed -abam tophat_out/accepted_hits_nojunc.bam -b ${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes_exons.bed -s -f 1.0 -wa -wb -bed > tophat_out/overlaps_exons.bed
    echo "intersectBed returned $?" >&2

    # Eliminate reads that map to more than a single gene
    awk '{ print $_"\t"$4; }' <tophat_out/overlaps_exons.bed | sort -k 13 | uniq -f 12 -u | cut -f 1-12 >tophat_out/overlaps_exons_uniq.bed

fi  # END !STATSONLY


# collect statistics
stats=($(awk '$2 ~ /^[0-9]+$/ {print $2 }' <clipstats.txt))
reads[raw]=${stats[0]}
reads[tooshort]=$(( ${reads[raw]} - ${stats[2]} ))
reads[adapter]=$(( ${reads[tooshort]} - ${stats[3]} ))
reads[uncalled]=$(( ${reads[adapter]} - ${stats[4]} ))


log_msg "Counting rRNA alignments..."
reads[rrna]=$(readc non-rrna.fq)

# Count how many failed to align
log_msg "Counting failed alignments..."
reads[unaligned]=$(( $(samtools view tophat_out/accepted_hits.bam | awk '{print $1}' |sort|uniq|wc -l) ))
alignments[unaligned]=$(samtools view tophat_out/accepted_hits.bam | wc -l )


# Count how many imperfect alignments
log_msg "Counting imperfect alignments..."
reads[imperfect]=$(($(samtools view tophat_out/accepted_hits_perfect.bam | awk '{print $1}' | sort|uniq|wc -l) ))

alignments[imperfect]=$(samtools view tophat_out/accepted_hits_perfect.bam | wc -l)

# Count how many alignments without junctions
log_msg "Counting alignments within a single exon..."
reads[nojunc]=$(($(samtools view tophat_out/accepted_hits_nojunc.bam | awk '{print $1}' | sort|uniq|wc -l) ))

alignments[nojunc]=$(samtools view tophat_out/accepted_hits_nojunc.bam | wc -l)

# how many aligned reads fell within known exons
log_msg "Counting reads that are contained within KNOWN exons..."
reads[exons]=$(awk '{print $4}' <tophat_out/overlaps_exons.bed | sort|uniq|wc -l)
alignments[exons]=$(wc -l <tophat_out/overlaps_exons.bed)

# exclude reads that overlap multiple features (genes) at a locus
log_msg "Counting overlaps with single features..."
reads[mult_gene]=$(awk '{ print $4; }' <tophat_out/overlaps_exons_uniq.bed | sort | uniq -c | wc -l )

alignments[mult_gene]=$(wc -l <tophat_out/overlaps_exons_uniq.bed)
    
#
# THIS IS IT!  This next call to rpos_dist.py is the heart of this
# routine and the main reason we have prepared the reads up to this point.
#
# Now that we know which gene each read maps to, and each mapping is
# unique, we can calculate where the read falls along the gene.  Not
# only do we get information about which region of the gene the read
# is in (5'-UTR, CDS, or 3'-UTR), we also get information about what
# reading frame the read is aligned with.  Since each read represents
# a ribosome-protected tag, we can derive information about ribosome
# density on genes.

log_msg "preparing coverage statistics..."
python ${SCRIPTS}/rpos_dist.py ${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes.txt tophat_out/overlaps_exons_uniq.bed  >tophat_out/aligned_position_stats.txt 2>rpos_dist_stats.txt
echo "python rpos_dist.py returned $?" >&2

log_msg "graphing results..."
${SCRIPTS}/plot_rpos_vs_region.sh <tophat_out/aligned_position_stats.txt
log_msg " "
sleep 3
echo 

cat flow.cmd
printf "%20s\t%d\t%d\n" "Raw reads" ${reads[raw]} ${alignments[raw]}
printf "%20s\t%d\t%d\n" "too short (<20nt)" ${reads[tooshort]} ${alignments[tooshort]}
printf "%20s\t%d\t%d\n" "adapter only" ${reads[adapter]} ${alignments[adapter]}
printf "%20s\t%d\t%d\n" "uncalled base 'N'" ${reads[uncalled]} ${alignments[uncalled]}
printf "%20s\t%d\t%d\n" "non-rRNA/non-spike" ${reads[rrna]} ${alignments[rrna]}
printf "%20s\t%d\t%d\n" "aligned" ${reads[unaligned]} ${alignments[unaligned]}
printf "%20s\t%d\t%d\n" "perfect alignment" ${reads[imperfect]} ${alignments[imperfect]}
printf "%20s\t%d\t%d\n" "nojunc alignment" ${reads[nojunc]} ${alignments[nojunc]}
printf "%20s\t%d\t%d\n" "mapped exons" ${reads[exons]} ${alignments[exons]}
printf "%20s\t%d\t%d\n" "single genes" ${reads[mult_gene]} ${alignments[mult_gene]}

echo
cat rpos_dist_stats.txt