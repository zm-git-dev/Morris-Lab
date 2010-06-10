#!/bin/bash

reads="$1"
outfile=${reads%.*}

bowtie -S -f mm_fmr1 "${reads}" "${outfile}".sam
samtools view -bS -o "${outfile}.bam" "${outfile}.sam"
samtools sort "${outfile}.bam" "${outfile}"
samtools index "${outfile}.bam"
echo "Produced \"${outfile}.bam\""
