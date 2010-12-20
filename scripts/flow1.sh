#!/bin/bash
# Workflow for aligning reads to the mm_fmr genome
#
# Usage:  flow1.sh <reads.fa>  
#

reads="$1"
outfile=${reads%.*}

bowtie -S --best -n 1 -f mm_fmr1 "${reads}" "${outfile}".sam
samtools view -bS -o "${outfile}.bam" "${outfile}.sam"
samtools sort "${outfile}.bam" "${outfile}"
samtools index "${outfile}.bam"
echo "Produced \"${outfile}.bam\""
