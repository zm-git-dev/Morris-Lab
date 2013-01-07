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

all: ${aligner}_out/aligned_position_stats.txt

clean: clean-common clean-bowtie clean-tophat

clean-common:
	-rm -f reads_trimmed.fq
	-rm -f rrna_reads.fq reads_nonrrna.fq
	-rm -f cutadapt.stats clipstats.txt
	-rm -f rrna.stats tooshort.fq rrna_aligned.sam rrna.bam

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
ifeq ($(REVERSECOMPLEMENT),1)
	${cmd} $^ | cutadapt -f fastq -m ${LENGTH} -a ${ADAPTER} --too-short-output=tooshort.fq /dev/stdin 2> cutadapt.stats | fastx_reverse_complement -Q33 -o $@
else
	${cmd} $^ | cutadapt -f fastq -m ${LENGTH} -a ${ADAPTER} --too-short-output=tooshort.fq /dev/stdin 2>&1 >$@ | tee cutadapt.stats
endif

$(BT)/bowtie_viral_hits.sam : $(BT)/bowtie_nomatch.fastq
	bowtie --best -k 3 -v ${MISMATCHES} --sam "${REFDIR}/${GENOME}/viruses"  $^ >$@

$(BT)/bowtie_nomatch.fastq : $(BT)/bowtie_hits.sam

#
# Align the reads.
# We can use either bowtie or tophat, depending on the target.
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

#
# filter out alignments that have more than ${MISMATCHES} base mismatches.
# FIXME: this fails if there are 10 or more mismatches
#
${aligner}_out/accepted_hits_perfect.sam: ${aligner}_out/accepted_hits.bam
	samtools view -h $^ | awk '/^@/ { print;} /NM:i:[0-9]+/ { if (int(substr($$0,match($$0, /NM:i:/)+5)) <= ${MISMATCHES}) print; }' >$@

#
# filter out reads that align to ribosomal DNA
#
reads_nonrrna.fq : reads_trimmed.fq
	bowtie -p 4 -S -n 2 -e 70 -l 28 --maxbts 125 -k 1 --al reads_rrna.fq --un $@ --phred33-quals "${REFDIR}/${GENOME}/rrna" $^  >rrna_aligned.sam 2>rrna.stats  

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
	intersectBed -abam $^ -b ${knowngenes}_exons.bed -s -f 1.0 -wa -wb -bed |cut -f 1-6,13- > $@

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
dbstats.txt:  ${aligner}_out/accepted_hits.bam ${aligner}_out/overlaps_exons.bed 
	@echo "making dbstats.txt..."
	@grep "Processed reads" cutadapt.stats | awk '{print $$3}' >$@
	@wc -l reads_trimmed.fq | awk '{print $$1/4}' >>$@
	@grep "failed to align" rrna.stats | awk '{print $$7}' >>$@
	@samtools view  ${aligner}_out/accepted_hits.bam | cut -f 1 | sort | uniq | wc -l  >>$@
#	@wc -l ${aligner}_out/bowtie_aligned.fastq | awk '{print $$1/4}' >>$@
	@cut -f 4 ${aligner}_out/overlaps_exons.bed | sed 's!/[0-9]!!' | sort | uniq | wc -l >>$@
	@cut -f 4 ${aligner}_out/overlaps_exons_uniq.bed | sed 's!/[0-9]!!' | sort | uniq | wc -l >>$@

db_stats:  dbstats.txt
	$(SCRIPTS)/sql_enter_stats -g "$(GENOME)" -m "$(MISMATCHES)" -d "$(DATASET)" "$(EXPERIMENT)" <dbstats.txt

db_expressions: ${aligner}_out/fpkm.out
	$(SCRIPTS)/sql_enter_expressions "$(DATASET)" "${aligner}_out/fpkm.out"

db_alignments:  ${aligner}_out/accepted_hits.bam 
	$(SCRIPTS)/sql_enter_alignments "$(DATASET)" "${aligner}_out/accepted_hits.bam"

database: db_stats db_alignments db_expressions

${aligner}_out/final_alignments.bam: ${aligner}_out/overlaps_exons_uniq.bed
	(samtools view -H ${aligner}_out/accepted_hits.bam; \
	 samtools view ${aligner}_out/accepted_hits.bam | fgrep -f <( cut -f 4 $< )) | \
		samtools view -S -b /dev/stdin > ${aligner}_out/igv.bam
	samtools sort ${aligner}_out/igv.bam ${aligner}_out/igv_s
	samtools index ${aligner}_out/igv_s.bam
	mv ${aligner}_out/igv_s.bam  ${aligner}_out/final_alignments.bam
	mv ${aligner}_out/igv_s.bam.bai  ${aligner}_out/final_alignments.bam.bai

rrna.bam : rrna_aligned.sam
	samtools view -S -b $< >$@
	samtools sort rrna.bam rrna_s
	samtools index rrna_s.bam
	mv rrna_s.bam rrna.bam
	mv rrna_s.bam.bai rrna.bam.bai

sample.bam : rrna.bam
	samtools view -b -s 0.001 rrna.bam >sample.bam
	samtools sort sample.bam sample_s
	samtools index sample_s.bam
	mv sample_s.bam sample.bam
	mv sample_s.bam.bai sample.bam.bai

HTML=/var/www/html
export EXPERIMENT

install-alignments: rrna.bam ${aligner}_out/final_alignments.bam
	-mkdir -p  $(HTML)/alignments/$${EXPERIMENT%_*}/${EXPERIMENT}
	cp ${aligner}_out/final_alignments.bam $(HTML)/alignments/$${EXPERIMENT%_*}/${EXPERIMENT}/${EXPERIMENT}.bam
	cp ${aligner}_out/final_alignments.bam.bai $(HTML)/alignments/$${EXPERIMENT%_*}/${EXPERIMENT}/${EXPERIMENT}.bam.bai
	cp rrna.bam $(HTML)/alignments/$${EXPERIMENT%_*}/${EXPERIMENT}/rrna_${EXPERIMENT}.bam
	cp rrna.bam.bai $(HTML)/alignments/$${EXPERIMENT%_*}/${EXPERIMENT}/rrna_${EXPERIMENT}.bam.bai

rrna_depth.txt: rrna.bam
	samtools depth rrna.bam | sort -k 3rn >rrna_depth.txt

.PHONY: clean clean-tophat clean-bowtie clobber clobber-tophat clobber-bowtie stats # dbstats.txt

.DELETE_ON_ERROR :
