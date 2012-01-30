
MISMATCHES=0
define perfect_filter
/^@/ {print}  
!/^@/ { for (i = 1; i <= NF; i++) 
	if ($$i ~ /\<NM:i:[0-9]+\>/ && match($$i, /[0-9]+/) && substr($$i, RSTART, RLENGTH) >= "${MISMATCHES}") { 
		print ; 
		break; 
	} 
}
endef

export perfect_filter

GENOME=mm9
INDEX_BASE=${GENOME}/${GENOME}
RRNA_BASE=${GENOME}/rrna
REFSEQ_BASE=${GENOME}
MORRIS=/home/csw/Morris-Lab
REFDIR=/usr/local/share/genome
SCRIPTS=${MORRIS}/scripts
KNOWNGENES=${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes

all: tophat_out/aligned_position_stats.txt

clean:

clobber:
	-rm -f tophat_out/aligned_position_stats.txt  tophat_out/overlaps_exons_uniq.bed tophat_out/accepted_hits_nojunc.bam tophat_out/accepted_hits_perfect.bam

plot: tophat_out/aligned_position_stats.txt
	${SCRIPTS}/plot_rpos_vs_region.sh <$^

tophat_out/accepted_hits.bam: bowtie_hits.sam
	-mkdir $(@D)
	awk '$$3!="*"' <$^ | samtools view -Sbh /dev/stdin >$@

tophat_out/accepted_hits_perfect.sam: tophat_out/accepted_hits.bam
	#samtools view -h $^ | grep -w 'NM:i:${MISMATCHES}' >$@
	samtools view -h $^  >$@

tophat_out/accepted_hits_nojunc.bam:  tophat_out/accepted_hits_perfect.sam
	awk '/^@/ || $$6 ~ /^[0-9]+M$$/' $^ | samtools view -Sbh /dev/stdin >$@

tophat_out/overlaps_exons.bed: tophat_out/accepted_hits_nojunc.bam
	intersectBed -abam $^ -b ${KNOWNGENES}_exons.bed -s -f 1.0 -wa -wb -bed > $@

tophat_out/overlaps_exons_uniq.bed: tophat_out/overlaps_exons.bed
	awk '{ print $$_"\t"$$4; }' <$^ | sort -k 13 | uniq -f 12 -u | cut -f 1-12 >$@

tophat_out/aligned_position_stats.txt: tophat_out/overlaps_exons_uniq.bed
	python ${SCRIPTS}/rpos_dist.py ${KNOWNGENES}.txt $^  >$@ 2>rpos_dist_stats.txt


stats: tophat_out/accepted_hits_perfect.sam tophat_out/accepted_hits_nojunc.bam tophat_out/overlaps_exons.bed tophat_out/overlaps_exons_uniq.bed
	@echo -n "perfect alignments:\t"
	#@echo -n $$(awk '!/^@/ {print $$1}' tophat_out/accepted_hits_perfect.sam | sort|uniq|wc -l)"\t"	
	@awk '!/^@/ {print $$1}' tophat_out/accepted_hits_perfect.sam | sort|uniq -c| awk '{ total += $$1 } END { print "perfect alignments:",NR,total; }'

	@echo -n "w/o junctions:\t"
	@echo -n $$(samtools view tophat_out/accepted_hits_nojunc.bam | awk '{print $$1}' | sort|uniq|wc -l)"\t"
	@echo $$(samtools view tophat_out/accepted_hits_nojunc.bam | wc -l)

	@echo -n "in exon:\t"
	@echo -n $$(awk '{print $$4}' <tophat_out/overlaps_exons.bed | sort|uniq|wc -l)"\t"
	@echo $$(wc -l <tophat_out/overlaps_exons.bed)

	@echo -n "one gene @ locus:\t"
	@echo -n $$(awk '{ print $$4; }' <  tophat_out/overlaps_exons_uniq.bed | sort | uniq -c | wc -l )"\t"
	@echo $$(wc -l <  tophat_out/overlaps_exons_uniq.bed)
