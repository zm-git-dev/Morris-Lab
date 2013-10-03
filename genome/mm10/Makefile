# Build genomes needed for alignment and analysis
#
# This makefile was written using information supplied by
# http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
SHELL=/bin/bash
FASTA_DIR=fa
GENOME=mm10

knowngenes: ${GENOME}_refseq_knowngenes.gtf

# copy the knownegenes table directly from UCSC database server.
knowngenes_db: ${GENOME}_refseq_knowngenes.txt refGene.sql
	sed 's/refGene/${GENOME}_knowngenes/g' <refGene.sql | mysql

${GENOME}_refseq_knowngenes.gtf: genePredToGtf  ${GENOME}_refseq_knowngenes.txt
	cut -f2-  ${GENOME}_refseq_knowngenes.txt | ./genePredToGtf file stdin ${GENOME}_refseq_knowngenes.gtf

# The genePredToGtf tool processes files in GenePred format into GTF format.
# See link to ucsc above for more information.   A version of this tools is 
# also available for Macs.
genePredToGtf:
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
	chmod +x genePredToGtf

${GENOME}_refseq_knowngenes.txt: 
	wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/database/refGene.txt.gz -O ${GENOME}_refseq_knowngenes.txt.gz
	gunzip -f  ${GENOME}_refseq_knowngenes.txt.gz

db_knowngenes: refGene.sql
	./sql_update_knowngenes ${GENOME} refGene.sql

refGene.sql : 
	mysqldump -h genome-mysql.cse.ucsc.edu -u genomep -ppassword --skip-lock-tables --tables ${GENOME} refGene >$@

# Build a genome index suitable for use with Bowtie and TopHat.
#
build: genome
	@echo "The build process can take a LONG, LONG TIME!"
	@echo "(as in several hours!)"
	@-read -t 10 -n 1 -p "Are you sure you want to continue? (Y/n) "
	bowtie-build $$(ls -m ${FASTA_DIR}/*.fa |  tr '\n' ' ' | sed 's/,\s*/,/g') ${GENOME}

# Retrieve the FASTA files used to build a genome.
#
genome:
	-mkdir ${FASTA_DIR}
	cd ${FASTA_DIR} && wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/chromosomes/chr[1-9XY].fa.gz'
	cd ${FASTA_DIR} && wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/chromosomes/chr[1-9]?.fa.gz'
	cd ${FASTA_DIR} && gunzip -f *.fa.gz

