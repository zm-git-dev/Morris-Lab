# Build genomes needed for alignment and analysis
#
# This makefile was written using information supplied by
# http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
SHELL=/bin/bash
SCRIPTS=~morrislab/scripts
FASTA_DIR=fa
FASTA_GZ=fa_gz
GENOME=hg18
chr_no   := 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
chr_id   := $(foreach x,$(chr_no),chr${x})
chr_fa   := $(foreach x,$(chr_no),chr${x}.fa)
chr_fai  := $(foreach x,$(chr_no),chr${x}.fa.fai)
empty:=
space:= $(empty) $(empty)
comma:=,
chr_pat   := $(subst $(space),|,$(foreach x,$(chr_no),chr${x}))
HOSTNAME := mauritius.bchem.washington.edu
WWW_GENOMES := igv/genomes
GENOMES_ROOT := /var/www/html/${WWW_GENOMES}
GENOME_DIR := ${GENOMES_ROOT}/${GENOME}

include ../genome.mk
