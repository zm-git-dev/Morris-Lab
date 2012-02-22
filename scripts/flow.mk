
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
ADAPTER=TGGAATTCTCGGGTGCCAAGG
LENGTH=20
REVERSECOMPLEMENT=0

all: ${aligner}_out/aligned_position_stats.txt

clean: clean-common clean-bowtie clean-tophat

clean-common:
	-rm -f ${basename}.multilen.fq
	-rm -f ${basename}.nonrrna.fq
	-rm -f cutadap.out clipstats.txt

clean-tophat: 
	-rm -f tophat_out/aligned_position_stats.txt  tophat_out/overlaps_exons_uniq.bed tophat_out/accepted_hits_nojunc.bam tophat_out/accepted_hits_perfect.bam

clean-bowtie: 
	-rm -f bowtie_out/aligned_position_stats.txt  bowtie_out/overlaps_exons_uniq.bed bowtie_out/accepted_hits_nojunc.bam bowtie_out/accepted_hits_perfect.bam

clobber: clobber-common clobber-bowtie clobber-tophat

clobber-tophat:
	-rm -rf tophat_out

clobber-bowtie: clobber-common
	-rm -rf bowtie_out

clobber-common : clean-common

plot: ${aligner}_out/aligned_position_stats.txt
	bash ${SCRIPTS}/plot_rpos_vs_region.sh <$^

cmd.gz=gunzip -c
cmd=$(if $(cmd$(suffix ${rawreads})),$(cmd$(suffix ${rawreads})),cat)

${basename}.multilen.fq : ${rawreads}
	#${cmd} $^ | fastx_clipper -a ${ADAPTER} -l ${LENGTH} -Q33 -v >$@ 
ifeq ($(REVERSECOMPLEMENT),1)
	${cmd} $^ | cutadapt -f fastq -m ${LENGTH} -a ${ADAPTER} /dev/stdin | fastx_reverse_complement -Q33 -o $@
else
	${cmd} $^ | cutadapt -f fastq -m ${LENGTH} -a ${ADAPTER} /dev/stdin 2>&1 >$@ | tee cutadapt.stats
endif

#
# Align the reads.
# We can use either bowtie or tophat, depending on the target.
#
bowtie_out/bowtie_hits.sam : ${basename}.nonrrna.fq
	@-mkdir -p $(dir $@)
	bowtie --best -k 3 -v ${MISMATCHES} --sam "/${INDEX_BASE}"  $^ >$@

bowtie_out/accepted_hits.bam: bowtie_out/bowtie_hits.sam
	#awk '/^@/ || $$3!="*"' <$^ | samtools view -Sbh /dev/stdin >$@
	${SCRIPTS}/sam_renamer <$^ | samtools view -Sbh /dev/stdin >$@

tophat_out/accepted_hits.bam: ${basename}.nonrrna.fq 
	@-mkdir -p $(dir $@)
	tophat --segment-length 10 --segment-mismatches ${MISMATCHES} -G  ${knowngenes}.gtf  "/${INDEX_BASE}" $^

#
# filter out alignments that have more than ${MISMATCHES} base mismatches.
# FIXME: this fails if there are 10 or more mismatches
#
${aligner}_out/accepted_hits_perfect.sam: ${aligner}_out/accepted_hits.bam
	samtools view -h $^ | awk '/^@/ { print;} /NM:i:[0-9]+/ { if (int(substr($$0,match($$0, /NM:i:/)+5)) <= ${MISMATCHES}) print; }' >$@

#
# filter out reads that align to ribosomal DNA
#
${basename}.nonrrna.fq : ${basename}.multilen.fq
	bowtie --un $@ "${RRNA_BASE}" $^ >/dev/null  


#
# filter out alignments that are not continuous (they span an intron).
# This really only applies to alignments from tophat for now, although
# in the future, bowtie II will also find gapped alignments.
#
# The 6th field in a SAM record is a CIGAR string that indicates how
# the read was aligned to the reference genome.  CIGAR strings of ##M
# (e.g. 21M or 36M, etc) indicates how many bases continuously matched
# or mismatched the genome.  A CIGAR string like 10M2003N15M indicates
# that there was a long string of unmatched bases in the middle
# ("2003N"), and that usually means there is a junction or intron in
# the gap.  We want to filter those out of our alignment files.
#
${aligner}_out/accepted_hits_nojunc.bam:  ${aligner}_out/accepted_hits_perfect.sam
	awk '/^@/ || $$6 !~ /.*[0-9]*M[0-9]+N[0-9]+M$$/' $^ | samtools view -Sbh /dev/stdin >$@

${aligner}_out/overlaps_exons.bed: ${aligner}_out/accepted_hits_nojunc.bam
	intersectBed -abam $^ -b ${knowngenes}_exons.bed -s -f 1.0 -wa -wb -bed > $@

${aligner}_out/overlaps_exons_uniq.bed: ${aligner}_out/overlaps_exons.bed
	awk '{ print $$_"\t"$$4; }' <$^ | sort -k 13 | uniq -f 12 -u | cut -f 1-12 >$@

${aligner}_out/aligned_position_stats.txt: ${aligner}_out/overlaps_exons_uniq.bed
	python ${SCRIPTS}/rpos_dist.py ${knowngenes}.txt $^ 2>&1  >$@ | tee ${aligner}_out/rpos_dist_stats.txt

#
# Targets for detecting ligase-bias
#
ligase-bias : ${aligner}_out/bestseq.txt
	echo "5' distribition"
	sed 's/\(.\).*/\1/' <$^ | sort | uniq -c | sort -k 2
	echo "3' distribution"
	sed 's/.*\(.\)/\1/' <$^ | sort | uniq -c | sort -k 2

${aligner}_out/bestseq.txt : ${aligner}_out/bestreads.txt  ${aligner}_out/accepted_hits.bam
	samtools view ${aligner}_out/accepted_hits.bam | fgrep -f  ${aligner}_out/bestreads.txt | cut -f 10 | grep -v N >$@

${aligner}_out/bestreads.txt : ${aligner}_out/aligned_position_stats.txt
	cut -f 1 $^ | sed 's/\/[0-9]*/\/1/' | head -50000 | sort | uniq >$@

#
# pseudo target to calculate various statistics about the workflow.
#
stats: read_stats alignment_stats

alignment_stats : ${aligner}_out/accepted_hits_perfect.sam ${aligner}_out/accepted_hits_nojunc.bam ${aligner}_out/overlaps_exons.bed ${aligner}_out/overlaps_exons_uniq.bed
	@echo -n "perfect alignments:\t"
	@echo -n $$(awk '!/^@/ {print $$1}' ${aligner}_out/accepted_hits_perfect.sam | sort|uniq|wc -l)"\t"	
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

.DELETE_ON_ERROR :
