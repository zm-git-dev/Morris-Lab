
rawreads=reads.fq
basename=$(basename $(notdir ${rawreads}))
MISMATCHES=0
aligner=tophat
GENOME=mm9
INDEX_BASE=${GENOME}/${GENOME}
RRNA_BASE=${GENOME}/rrna
REFSEQ_BASE=${GENOME}
MORRIS=/home/csw/Morris-Lab
REFDIR=/usr/local/share/genome
SCRIPTS=${MORRIS}/scripts
knowngenes=${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes
ADAPTER="ATCTCGTATGCCGTCTTCTGCTTG"
LENGTH=20

all: ${aligner}_out/aligned_position_stats.txt

clean: clean-bowtie clean-tophat

clean-tophat: 
	-rm -f tophat_out/aligned_position_stats.txt  tophat_out/overlaps_exons_uniq.bed tophat_out/accepted_hits_nojunc.bam tophat_out/accepted_hits_perfect.bam

clean-bowtie: 
	-rm -f bowtie_out/aligned_position_stats.txt  bowtie_out/overlaps_exons_uniq.bed bowtie_out/accepted_hits_nojunc.bam bowtie_out/accepted_hits_perfect.bam

clobber: clobber-bowtie clobber-tophat

clobber-tophat:
	-rm -f tophat_out/aligned_position_stats.txt  tophat_out/overlaps_exons_uniq.bed tophat_out/accepted_hits_nojunc.bam tophat_out/accepted_hits_perfect.bam

clobber-bowtie: clobber-common
	-rm -f bowtie_out/bowtie_hits.sam bowtie_out/aligned_position_stats.txt  bowtie_out/overlaps_exons_uniq.bed bowtie_out/accepted_hits_nojunc.bam bowtie_out/accepted_hits_perfect.bam

clobber-common : 
	-rm -f short.multilen.fastq

plot: ${aligner}_out/aligned_position_stats.txt
	${SCRIPTS}/plot_rpos_vs_region.sh <$^

${basename}.multilen.fastq : ${rawreads}
	fastx_clipper -a "${ADAPTER}" -l ${LENGTH} -Q33 -v <$^ > $@ 

#
# Align the reads.
# We can use either bowtie or tophat, depending on the target.
#
bowtie_out/bowtie_hits.sam : bowtie_out/non-rrna.fq 
	bowtie --best --sam "/${INDEX_BASE}"  $^ >$@

bowtie_out/accepted_hits.bam: bowtie_out/bowtie_hits.sam
	awk '/^@/ || $$3!="*"' <$^ | samtools view -Sbh /dev/stdin >$@

tophat_out/accepted_hits.bam: reads.multilen.fastq
	tophat --segment-length 10 --segment-mismatches 0 -G  ${knowngenes}.gtf  "/${INDEX_BASE}" $^

#
# filter out alignments that have more than ${MISMATCHES} base mismatches.
# FIXME: this fails if there are 10 or more mismatches
#
${aligner}_out/accepted_hits_perfect.sam: ${aligner}_out/accepted_hits.bam
	samtools view -h $^ | grep -w 'NM:i:[0-${MISMATCHES}]' >$@

#
# filter out reads that align to ribosomal DNA
#
${aligner}_out/non-rrna.fq : ${basename}.multilen.fastq
	@-mkdir $(@D)
	bowtie --un $@ "${RRNA_BASE}" ${basename}.multilen.fastq >/dev/null  


#
# filter out alignments that are not continuous (they span an intron).
# This really only applies to alignments from tophat for now, although
# in the future, bowtie II will also find gapped alignments.
#
${aligner}_out/accepted_hits_nojunc.bam:  ${aligner}_out/accepted_hits_perfect.sam
	awk '/^@/ || $$6 ~ /^[0-9]+M$$/' $^ | samtools view -Sbh /dev/stdin >$@

${aligner}_out/overlaps_exons.bed: ${aligner}_out/accepted_hits_nojunc.bam
	intersectBed -abam $^ -b ${knowngenes}_exons.bed -s -f 1.0 -wa -wb -bed > $@

${aligner}_out/overlaps_exons_uniq.bed: ${aligner}_out/overlaps_exons.bed
	awk '{ print $$_"\t"$$4; }' <$^ | sort -k 13 | uniq -f 12 -u | cut -f 1-12 >$@

${aligner}_out/aligned_position_stats.txt: ${aligner}_out/overlaps_exons_uniq.bed
	python ${SCRIPTS}/rpos_dist.py ${knowngenes}.txt $^  >$@ 2>rpos_dist_stats.txt

#
# Calculate various statistics about the workflow.
#
stats: ${aligner}_out/accepted_hits_perfect.sam ${aligner}_out/accepted_hits_nojunc.bam ${aligner}_out/overlaps_exons.bed ${aligner}_out/overlaps_exons_uniq.bed
	@echo -n "perfect alignments:\t"
	#@echo -n $$(awk '!/^@/ {print $$1}' ${aligner}_out/accepted_hits_perfect.sam | sort|uniq|wc -l)"\t"	
	@awk '!/^@/ {print $$1}' ${aligner}_out/accepted_hits_perfect.sam | sort|uniq -c| awk '{ total += $$1 } END { print "perfect alignments:",NR,total; }'

	@echo -n "w/o junctions:\t"
	@echo -n $$(samtools view ${aligner}_out/accepted_hits_nojunc.bam | awk '{print $$1}' | sort|uniq|wc -l)"\t"
	@echo $$(samtools view ${aligner}_out/accepted_hits_nojunc.bam | wc -l)

	@echo -n "in exon:\t"
	@echo -n $$(awk '{print $$4}' <${aligner}_out/overlaps_exons.bed | sort|uniq|wc -l)"\t"
	@echo $$(wc -l <${aligner}_out/overlaps_exons.bed)

	@echo -n "one gene @ locus:\t"
	@echo -n $$(awk '{ print $$4; }' <  ${aligner}_out/overlaps_exons_uniq.bed | sort | uniq -c | wc -l )"\t"
	@echo $$(wc -l <  ${aligner}_out/overlaps_exons_uniq.bed)

.PHONY: clean clean-tophat clean-bowtie clobber clobber-tophat clobber-bowtie stats

