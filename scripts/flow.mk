SHELL := /bin/bash
rawreads ?= reads.fq
MISMATCHES ?= 0
aligner ?= bowtie
GENOME ?= mm9
ADAPTER ?= TGGAATTCTCGGGTGCCAAGG
LENGTH ?= 20
REVERSECOMPLEMENT ?= 0
BOWTIE_INDEXES ?= /home/morrislab/genome
export BOWTIE_INDEXES

basename=$(basename $(notdir ${rawreads}))
INDEX_BASE=${GENOME}/${GENOME}
RRNA_BASE=${GENOME}/rrna
REFSEQ_BASE=${GENOME}
MORRIS ?= /home/morrislab
REFDIR ?= /home/morrislab/genome
SCRIPTS=${MORRIS}/scripts
knowngenes=${REFDIR}/${REFSEQ_BASE}/${REFSEQ_BASE}_refseq_knowngenes
TH=tophat_out
BT=bowtie_out

all: ${aligner}_out/accepted_final.sam

clean: clean-common clean-bowtie clean-tophat

clean-common:
	-rm -f reads_trimmed.fq
	-rm -f reads_rrna.fq reads_nonrrna.fq
	-rm -f cutadapt.stats clipstats.txt
	-rm -f rrna.stats tooshort.fq rrna_aligned.sam 
	-rm -f rrna.bam rrna.bam.bai rrna_depth.txt
	-rm -f dbstats.txt
	-rm -f sample.bam sample.bam.bai

clean-tophat: 
	-rm -f	$(TH)/aligned_position_stats.txt	\
		$(TH)/overlaps_exons_uniq.bed 		\
		$(TH)/accepted_hits_nojunc.bam 		\
		$(TH)/accepted_hits_perfect.bam

clean-bowtie: 
	-rm -f	$(BT)/aligned_position_stats.txt	\
		$(BT)/overlaps_exons_uniq.bed		\
		$(BT)/accepted_hits_nojunc.bam		\
		$(BT)/accepted_hits_perfect.bam $(BT)/igv.bam

clobber: clobber-common clobber-bowtie clobber-tophat

clobber-tophat:
	-rm -rf $(TH)

clobber-bowtie: clobber-common
	-rm -rf $(BT)

clobber-common : clean-common
	-rm -f rrna_s.bam rrna_s.bam.bai

plot: ${aligner}_out/aligned_position_stats.txt
	bash ${SCRIPTS}/plot_rpos_vs_region.sh <$^

cmd.gz=gunzip -c
cmd=$(if $(cmd$(suffix ${rawreads})),$(cmd$(suffix ${rawreads})),cat)

reads_trimmed.fq tooshort.fq  : ${rawreads}
	${cmd} $^ | fastq_quality_trimmer -Q 33 -t 25 -l ${LENGTH} -i /dev/stdin -o $(basename $^)_trim.fq
ifeq ($(REVERSECOMPLEMENT),1)
	cat $(basename $^)_trim.fq | cutadapt -f fastq -m ${LENGTH} -a ${ADAPTER} --too-short-output=tooshort.fq /dev/stdin 2> cutadapt.stats | fastx_reverse_complement -Q33 -o $@
else
	cat $(basename $^)_trim.fq | cutadapt -f fastq -m ${LENGTH} -a ${ADAPTER} --too-short-output=tooshort.fq /dev/stdin 2>&1 >$@ | tee cutadapt.stats
endif

$(BT)/bowtie_viral_hits.sam : $(BT)/bowtie_nomatch.fastq
	bowtie --best -k 3 -v ${MISMATCHES} --sam "${REFDIR}/${GENOME}/viruses"  $^ >$@

$(BT)/bowtie_nomatch.fastq : $(BT)/bowtie_hits.sam

#
# Align the reads.
# We can use either bowtie or tophat, depending on the target.
#
# The RPT reads are aligned against a reference genome.
# Those read that don't align the first time are trimmed and given a second chance at alignment.
#
$(BT)/bowtie_hits.sam : reads_nonrrna.fq
	@-mkdir -p $(dir $@)
	bowtie -p 4 --un $(BT)/bowtie_nomatch.fastq --al $(BT)/bowtie_aligned.fastq --best -k 3 -v ${MISMATCHES} --sam "${REFDIR}/${GENOME}/${GENOME}"  $^ 2>&1 >$@ | tee $(BT)/bowtie_hits.stats 


$(BT)/accepted_hits.bam: $(BT)/bowtie_hits.sam
	# ${SCRIPTS}/sam_renamer <$^ | samtools view -Sbh /dev/stdin >$@
#	awk '/^@/ { print; } !/^@/ && $$3!="*" { print | " sort -k 1 " }' $^ |samtools view -Sbh $^ >$@
	awk '/^@/ { print; } !/^@/ && $$3!="*" { print  }' $^ | samtools view -Sbh /dev/stdin >$@

$(TH)/accepted_hits.bam: reads_nonrrna.fq 
	@-mkdir -p $(dir $@)
	BOWTIE_INDEXES="${BOWTIE_INDEXES}" tophat --segment-length 10 --segment-mismatches ${MISMATCHES} -G  ${knowngenes}.gtf  "/${INDEX_BASE}" $^



# exclude alignments that are aligned to more than one gene, or to no gene at all.
# filter out alignments with more than ${MISMATCHES} mismatched bases.
# filter out alignments that are aligned to more than one locus on the genome.
# exclude alignments that are aligned to more than one gene, or to no gene at all.

${aligner}_out/accepted_final.sam : ${aligner}_out/accepted_hits.bam
	# filter out alignments with more than ${MISMATCHES} mismatched bases.
	# filter out alignments that are aligned to more than one locus on the genome
	samtools view $^ | awk '/NM:i:[0-9]+/ { if (int(substr($$0,match($$0, /NM:i:/)+5)) <= ${MISMATCHES}) print; }'  \
	                 | awk '{ print $$0"\t"$$1; }' | uniq -f 14 -u | cut -f 1-14  >${aligner}_out/tmp.sam
	#
	# find overlaps with features.
	#
	htseq-count --mode=intersection-strict -i transcript_id -o ${aligner}_out/overlaps_exons.sam  ${aligner}_out/tmp.sam \
			"${REFDIR}/${GENOME}/${GENOME}_refseq_knowngenes.gtf" \
	| awk '$$2 != 0 ' >${aligner}_out/gene_count.txt
	#
	# exclude alignments that are aligned to more than one gene, or to no gene at all.
	#
	(samtools view -H $^; awk '$$15 !~ "ambiguous|no_feature"'  ${aligner}_out/overlaps_exons.sam) >$@
	rm -f ${aligner}_out/tmp.sam # ${aligner}_out/tmp2.sam

#
# filter out reads that align to ribosomal DNA
#
reads_nonrrna.fq : reads_trimmed.fq
	bowtie -p 4 -S -n 2 -e 70 -l 28 --maxbts 125 -k 1 --al reads_rrna.fq --un $@ --phred33-quals "${REFDIR}/${GENOME}/rrna" $^  >rrna_aligned.sam 2>rrna.stats  

#
# filter out alignments that are not continuous (they span an intron).
# This really only applies to alignments from tophat for now, although
# in the future, (I think) bowtie II will also find gapped alignments.
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
	intersectBed -abam $^ -b ${knowngenes}_exons.bed -s -f 1.0 -wa -wb -bed |cut -f 1-6,13- > $@

# remove any alignment that overlaps the genome in more than one place
# (e.g. on chr19 and on chr3), -OR- which overlaps more than one
# feature at a single locus (e.g. at chr19:1002002 overlaps two
# isoforms of the same gene, or different genes.)
#
# This is a really tricky filter!
# Basically we are filtering any duplicate alignments for the same
# read.  This has the effect of filtering for aligments to multiple
# loci, and for multiple alignments to the same locus but different
# features (genes).  Its pretty tricky because a tool like
# htseq-count, while filtering for ambiguous alignments to the same
# locus, will not filter for multiple aligmnents to different loci.
# You'll need to use a different tool if you use htseq-count.
#
# Update: No, the above explanation is slightly wrong.  Only duplicate
# alignments that overlap multiple *features* (e.g. gene exons) are
# removed here.  If the original alignment was ambiguous (aligned to
# multiple loci in the genome), it will not be removed if only one of
# those loci overlaps a feature!  This is probably wrong and not what
# we want.


${aligner}_out/overlaps_exons_uniq.bed: ${aligner}_out/overlaps_exons.bed
	awk '{ print $$0"\t"$$4; }' <$^ | sort -k 13 | uniq -f 12 -u | cut -f 1-12 >$@

${aligner}_out/aligned_position_stats.txt: ${aligner}_out/overlaps_exons_uniq.bed
	python ${SCRIPTS}/rpos_dist.py ${knowngenes}.txt $^ 2>&1  >$@ | tee ${aligner}_out/rpos_dist_stats.txt

#
# Targets for detecting ligase-bias
#
raw-ligase-bias :  reads_trimmed.fq
	paste <($5primebias) <($3primebias)

	@echo "5' distribution"
	@${SCRIPTS}/reads <$^ |  grep -v N | sed 's/\(.\).*/\1/' | sort | uniq -c | sort -k 2
	@echo "3' distribution"
	@${SCRIPTS}/reads <$^ |  grep -v N | sed 's/.*\(.\)/\1/' | sort | uniq -c | sort -k 2

raw-rc-ligase-bias :  reads_trimmed.fq
	@echo "5' distribution"
	@${SCRIPTS}/reads <$^ |  grep -v N | ${SCRIPTS}/rc.pl | sed 's/\(.\).*/\1/' | sort | uniq -c | sort -k 2
	@echo "3' distribution"
	@${SCRIPTS}/reads <$^ |  grep -v N | ${SCRIPTS}/rc.pl | sed 's/.*\(.\)/\1/' | sort | uniq -c | sort -k 2




# Calculate a distribution of which remnants are too short.
#
# After the adapter has been stripped from the 3'-end of the reads,
# many reads are too short to process further (usually less than
# 20nt).  This rule calculates the variation among remaining
# fragments after stripping the adapter.  Fragments that are
# completely empty (length=0) represent adapter-dimers.
#
# Note: To speed calculation, only the first million fragments are considered.
# 

too-short : tooshort.fq
	# Look at the first million fragments after stripping adapter.
	$(SCRIPTS)/reads tooshort.fq | head -1000000 | sort | uniq -c | sort -rn 



# there is a subtle dependency here on bowtie_hits.sam, that will not
# work with accepted_hits.bam.  accepted_hit has been run through the
# renamer (See above) so each read name is unique.  The database entry
# goes to great lengths to correctly map the read names to a unique ID
# so you can see which alignments came from which reads.  Using
# accepted_hits will prevent this from working so just use
# bowtie_hits.sam for now.

 ${aligner}_out/fpkm.out :  ${aligner}_out/aligned_position_stats.txt
	$(SCRIPTS)/calculate_fpkm.py <$^ >$@


# prepare some statistics to be loaded into the database
# raw reads
# trimmed reads
# nonrrna reads
# aligned reads
# in exons
# in unique exons
dbstats.txt:  ${aligner}_out/accepted_hits.bam ${aligner}_out/accepted_final.sam
	@echo "making dbstats.txt..."
	@grep "Processed reads" cutadapt.stats | awk '{print $$3}' >$@
	@wc -l reads_trimmed.fq | awk '{print $$1/4}' >>$@
	@grep "failed to align" rrna.stats | awk '{print $$7}' >>$@
	@samtools view  ${aligner}_out/accepted_hits.bam | cut -f 1 | sort | uniq | wc -l  >>$@
#	@wc -l ${aligner}_out/bowtie_aligned.fastq | awk '{print $$1/4}' >>$@
	@grep -v '^@' ${aligner}_out/overlaps_exons.sam | wc -l >>$@
	@grep -v '^@' ${aligner}_out/accepted_final.sam | wc -l >>$@

db_stats:  dbstats.txt
	$(SCRIPTS)/sql_enter_stats -g "$(GENOME)" -m "$(MISMATCHES)" -d "$(DATASET)" "$(EXPERIMENT)" <dbstats.txt

db_alignments:  ${aligner}_out/accepted_final.sam 
	$(SCRIPTS)/sql_enter_new_alignments "$(DATASET)" "${aligner}_out/accepted_final.sam"

database: db_stats db_alignments

 ${aligner}_out/htseq_input.sam : ${aligner}_out/accepted_hits.bam
	samtools view $< | awk '!/^@/ { print $$0"\t"$$1; }' | uniq -f 14 -u | cut -f 1-14 >$@

 ${aligner}_out/htseq_hits.sam : ${aligner}_out/htseq_input.sam
	htseq-count --mode=intersection-strict -o $@ $< ~morrislab/genome/mm9/mm9_refseq_knowngenes.gtf 

${aligner}_out/htseq.bam: ${aligner}_out/htseq_hits.sam
	(samtools view -H ${aligner}_out/accepted_hits.bam; \
	 grep -v 'ambiguous\|no_feature' <$<) | \
		samtools view -S -b /dev/stdin > ${aligner}_out/htseq.bam
	samtools sort ${aligner}_out/htseq.bam ${aligner}_out/htseq_s
	samtools index ${aligner}_out/htseq_s.bam
	mv ${aligner}_out/htseq_s.bam  ${aligner}_out/htseq.bam
	mv ${aligner}_out/htseq_s.bam.bai  ${aligner}_out/htseq.bam.bai

igv: ${aligner}_out/accepted_final.bam ${aligner}_out/genomic.bam

${aligner}_out/accepted_final.bam: ${aligner}_out/accepted_final.sam
	samtools view -S -b  ${aligner}_out/accepted_final.sam > ${aligner}_out/igv.bam	
	samtools sort ${aligner}_out/igv.bam ${aligner}_out/igv_s
	samtools index ${aligner}_out/igv_s.bam
	mv ${aligner}_out/igv_s.bam  ${aligner}_out/accepted_final.bam
	mv ${aligner}_out/igv_s.bam.bai  ${aligner}_out/accepted_final.bam.bai

${aligner}_out/genomic.bam :  ${aligner}_out/accepted_hits.bam
	samtools sort ${aligner}_out/accepted_hits.bam ${aligner}_out/genomic_s
	samtools index ${aligner}_out/genomic_s.bam
	mv ${aligner}_out/genomic_s.bam ${aligner}_out/genomic.bam
	mv ${aligner}_out/genomic_s.bam.bai ${aligner}_out/genomic.bam.bai

# make a bam file of all the reads that aligned to ribosomal RNA so they can be
# displayed in IGV.  Note this may be too many reads to practically display in IGV.
# Use sample.bam to choose a random subsample of these reads.
rrna.bam : rrna_aligned.sam
	samtools view -S -b $< >$@
	samtools sort rrna.bam rrna_s
	samtools index rrna_s.bam
	mv rrna_s.bam rrna.bam
	mv rrna_s.bam.bai rrna.bam.bai

# Choose a random subsample (1/10000) of reads that aligned to
# ribosomal RNA so they can be displayed in IGV. This is the same data
# as in rrna.bam, but a smaller sample so IGV can easily handle it.
sample.bam : rrna.bam
	samtools view -b -s 0.001 rrna.bam >sample.bam
	samtools sort sample.bam sample_s
	samtools index sample_s.bam
	mv sample_s.bam sample.bam
	mv sample_s.bam.bai sample.bam.bai

HTML=/var/www/html
export EXPERIMENT

install-alignments: ${aligner}_out/genomic.bam ${aligner}_out/accepted_final.bam
	-mkdir -p  $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}
	-rm -f  $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/*.bam
	-rm -f  $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/*.bai
	cp ${aligner}_out/genomic.bam $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/${DATASET}_genomic.bam
	cp ${aligner}_out/genomic.bam.bai $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/${DATASET}_genomic.bam.bai
	cp ${aligner}_out/accepted_final.bam $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/${DATASET}_exons.bam
	cp ${aligner}_out/accepted_final.bam.bai $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/${DATASET}_exons.bam.bai
	# cp rrna.bam $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/${DATASET}_rrna.bam
	# cp rrna.bam.bai $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/${DATASET}_rrna.bam.bai
	# cp sample.bam $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/${DATASET}_smallrrna.bam
	# cp sample.bam.bai $(HTML)/alignments/$${EXPERIMENT%_*}/${DATASET}/${DATASET}_smallrrna.bam.bai

rrna_depth.txt: rrna.bam
	samtools depth rrna.bam | sort -k 3rn >rrna_depth.txt

.PHONY: clean clean-tophat clean-bowtie clobber clobber-tophat clobber-bowtie stats # dbstats.txt

.DELETE_ON_ERROR :
