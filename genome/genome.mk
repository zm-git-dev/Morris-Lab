# This makefile was written using information supplied by
# http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
SHELL=/bin/bash
SCRIPTS=~morrislab/scripts
FASTA_DIR=fa
FASTA_GZ=fa_gz
HOSTNAME := mauritius.bchem.washington.edu
WWW_GENOMES := igv/genomes
GENOMES_ROOT := /var/www/html/${WWW_GENOMES}
GENOME_DIR := ${GENOMES_ROOT}/${GENOME}

# define chr_no to include the chromosomes you wish you index.
# chr_no   := 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y

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

all: index ${GENOME}_local.genome

#
# makefile rules for building genome index files needed by Bowtie and IGV
#
# These rules must be specialized by defining make variables:
#

install: ${GENOMES_ROOT}
install: $(addprefix ${FASTA_DIR}/,$(chr_fai))  
install: ${GENOME}_local.genome
	-mkdir ${GENOME_DIR}
	tar cf - fa | (cd ${GENOME_DIR} && tar xf - )
	tar cf - ${GENOME}_local.genome | (cd ${GENOME_DIR} && tar xf - )

${GENOME}_local.genome: property.txt refGene.txt cytoBand.txt
	zip $@ $^

## The format of the property.txt file is defined by IGV
## this is the file that tell IGV where the genome sequence files are and 
## where the gene annotations can be found.
##
property.txt: Makefile
	echo "fasta=true" >$@
	echo "fastaDirectory=true" >>$@
	echo "fastaFiles=$(subst $(space),$(comma),${chr_fa})" >>$@
	echo "ordered=false" >>$@
	echo "id=${GENOME}_local" >>$@
	echo "name=${GENOME} local" >>$@
	echo "geneTrackName=RefSeq genes" >>$@
	echo "geneFile=refGene.txt" >>$@
	echo "cytobandFile=cytoBand.txt" >>$@
	echo "sequenceLocation=http://${HOSTNAME}/igv/genomes/${GENOME}/fa" >>$@
	echo "" >>$@


knowngenes: ${GENOME}_refseq_knowngenes.gtf

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

%.txt.gz :
	wget --quiet --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/database/$@

%.txt : %.txt.gz
	gunzip -c $^ >$@


# The refseq known genes may mention the same gene for multiple
# haplotypes, or it may reference genes on chromosomes that we are not
# interested in.   
# Filter out any mention of chromosomes that we do not want.
#
${GENOME}_refseq_knowngenes.txt: refGene.txt # $(GENOME)
	awk ' $$3 ~ /^(${chr_pat})$$/ { print }' <$^ >$@

# this rule is no longer used.
# I keep it around to remind myself how to fetch the table from the UCSC SQL 
# server in case I should want to do that again.
#
# refGene.sql : 
# 	mysqldump -h genome-mysql.cse.ucsc.edu -u genomep -ppassword --skip-lock-tables --tables ${GENOME} refGene >$@

# Retrieve the FASTA files used to build a genome.
#

genome: $(addprefix ${FASTA_DIR}/,$(chr_fai))


# Build a genome index suitable for use with Bowtie and TopHat.
#
index : $(GENOME).1.ebwt

$(GENOME).1.ebwt : $(addprefix ${FASTA_DIR}/,$(chr_fa))
	@echo "The build process can take a LONG, LONG TIME!"
	@echo "(as in several hours!)"
	bowtie-build $(subst $(space),$(comma),$^) $(GENOME)


${FASTA_DIR}/%.fa : ${FASTA_GZ}/%.fa.gz
	gunzip -c $<	>$@

.PRECIOUS: ${FASTA_GZ}/%.fa.gz

${FASTA_GZ}/%.fa.gz : ${FASTA_GZ} 
	cd ${FASTA_GZ} && wget --quiet --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/chromosomes/$(notdir $@)'

${FASTA_DIR} :
	mkdir -p ${FASTA_DIR} 

${FASTA_GZ} :
	mkdir -p ${FASTA_GZ} 


${FASTA_DIR}/%.fa.fai : ${FASTA_DIR}/%.fa
	samtools faidx $<

www:
	mkdir /var/www/genome/${GENOME}
	cp ${FASTA_DIR}/*.{fa,fai} /var/www/genome/${GENOME}/