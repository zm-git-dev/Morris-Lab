
# define chr_no to include the chromosomes you wish you index.
chr_no   := 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y

# prefix "chr" to the front of each chromosome identifier
chr_id   = $(foreach x,$(chr_no),chr${x})

# make FASTA fielnames from each chromosome identifier.
chr_fa   = $(foreach x,$(chr_no),chr${x}.fa)

# make FASTA index filenames from each chromosome identifier.
chr_fai  = $(foreach x,$(chr_no),chr${x}.fa.fai)

# make a comma-seperated list of chromosome identifiers
# (passed to bowtie-index on the command line.)
empty:=
space:= $(empty) $(empty)
comma:=,
chr_pat = $(subst $(space),|,$(foreach x,$(chr_no),chr${x}))

help:
	echo "make hg19"
	echo "make hg18"

hg19:
	make GENOME=hg19 hg19_union.gtf

hg18:
	make GENOME=hg18 hg18_union.gtf

mm9:
	make GENOME=mm9 mm9_union.gtf

mm10:
	make GENOME=mm10 mm10_union.gtf


${GENOME}_union.gtf: ${GENOME}_refseq_knowngenes.gtf
	Rscript union.R --genome "${GENOME}" $^

# copy the knownegenes table directly from UCSC database server.
knowngenes_db: ${GENOME}_refseq_knowngenes.txt
	${SCRIPTS}/sql_enter_knowngenes ${GENOME} $<

${GENOME}_refseq_knowngenes.gtf: genePredToGtf  ${GENOME}_refseq_knowngenes.txt
	cut -f2-  ${GENOME}_refseq_knowngenes.txt | ./genePredToGtf file stdin ${GENOME}_refseq_knowngenes.gtf

# The genePredToGtf tool processes files in GenePred format into GTF format.
# See link to ucsc above for more information.   A version of this tools is 
# also available for Macs.
genePredToGtf:
	wget --quiet http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
	chmod +x genePredToGtf

${GENOME}_%.txt.gz :
	wget --quiet --timestamping -O "${@}" ftp://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/database/$(patsubst ${GENOME}_%,%,${@})

%.txt : %.txt.gz
	gunzip -c $^ >$@


# The refseq known genes may mention the same gene for multiple
# haplotypes, or it may reference genes on chromosomes that we are not
# interested in.   
# Filter out any mention of chromosomes that we do not want.
#
${GENOME}_refseq_knowngenes.txt: ${GENOME}_refGene.txt # $(GENOME)
	awk ' $$3 ~ /^(${chr_pat})$$/ { print }' <$^ >$@


