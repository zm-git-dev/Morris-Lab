#!/usr/bin/env Rscript
##' Calculates pvalues for RNA-SEQ profiles from chemosensitive and chemoresistive phenotypes of ovarian tumor.
##'
##' @usage   chemotype [--help] [--gene <gene>] [--
##' @author  Chris Warth
##' @date 11-11-2013

suppressMessages(library(optparse))

# if executing locally, shared directories will be mounted under "/Volumes".
# if on the rhino servers, just go to the home directory.
prefix <- file.path("", "Volumes", "homes")
if (!file.exists(file.path(prefix, "data"))) {
    prefix <- file.path("~")
}

if (interactive()) {
    samplesize <- 1
    genes.file <- file.path(prefix, "data", "SOC", "genelist.txt")
    genes.additional <- c("WFDC2")
} else {
    option_list <- list(
        make_option("--normalization", default="area",
                    help = "Function to normalize data, \"quantile\" or \"area\" [default \"%default\"]"),
        make_option("--genefile", default="",
                    help="file containing name genes to monitor. [default \"%default\"]"),
        make_option("--sample", default="",
                    help="number of genes to sample from genefile. [default \"%default\"]"),
        make_option("--genes", default="",
                    help="gene or genes to monitor. [default \"%default\"]")
        )

    # get command line options, if help option encountered print help and exit,
    # otherwise if options not found on command line then set defaults,
    opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

    samplesize <- as.integer(opt$options$sample)
    if (is.na(samplesize))
        samplesize <- 0
    genes.file <- opt$options$genefile
    # split the list of genes, seperated by commas, and eliminate any empty strings
    genes.split <- gsub("\\s+", "", unlist(strsplit(opt$options$genes, ",")))
    genes.additional <- genes.split[nchar(genes.split)> 0]
}

suppressWarnings(suppressMessages(source("qarputils.R")))


## read SOC clinical outcomes spreadsheet
file.survival <- file.path(prefix, "data", "SOC", "2013-09-25_McIntoshSOC_Survival.csv")
df.survival <- read.csv(file.survival, colClasses=c(Sample="character"))
                      
## collect two cohorts of patient record that represent different treatment responses.
df.resistant <- subset(df.survival, Response == "Chemoresistant" & Resection == "Optimal")
df.sensitive <- subset(df.survival, Response == "Chemosensitive" & Resection == "Optimal")

cohort1 <- df.resistant[1:10,]
cohort2 <- df.sensitive[1:10,]

# Turn sample IDs into filenames
files <-  paste0(gsub("_","-", c(cohort1$Sample, cohort2$Sample)), ".bam")
files <- file.path(file.path(prefix, "data", "SOC"), files)

# The order of factors in aVec should mirror the order of BAM files
aVec <- as.factor(c(rep("chemoresistant", nrow(cohort1)), rep("chemosensitive", nrow(cohort2))))

# now read the GTF file to get descriptions of gene regions
gtfname <- system.file("python", "hg19-refSeq-union.gtf", package="qarp")
stopifnot(nchar(gtfname)>0)

gtf <- rtracklayer::import(gtfname, genome="hg19", feature.type=c("exon"), asRangedData=FALSE, format="gtf")
# delete non-canonical chromosomes (haplotypes and such)
# only keep features whose chromosome labels don't have underscores.
# e.g.  we throw out "chr3_hap6" but keep "chr3"
# Note: we should probably get rid of chrM, mitochondrial features, as well.
gtf <- gtf[seqnames(gtf) %in% levels(seqnames(gtf))[grep("_|chrM", levels(seqnames(gtf)), invert=TRUE)]]

set.seed(1234)
genes <- c()
if (nchar(genes.file) > 0) {
    genelist <- read.csv(genes.file, sep="\t", header=FALSE, stringsAsFactors = FALSE)
    if (samplesize > 0) {
        genes <- sample(genelist[, 1], samplesize)
    } else {
        genes <- genelist[,1]
    }
}

## split the additional genes list by comms followed by optoinal spaces
if (!all(nchar(genes.additional)==0))
    genes <- unique(c(genes.additional, genes))

totals <- sapply(files, fastBamCount)

# output the header.
# output the header.
cat(paste0("# BAM files: ", paste0(files, collapse=", "), "\n"))
cat(paste0("# phenotypes: ", paste0(aVec, collapse=", "), "\n"))
cat(paste0("# total alignments: ", paste0(totals, collapse=", "), "\n"))
cat(paste0("# annotations: ", gtfname, "\n"))
cat(paste0("# gene count: ", length(genes), "\n"))
cat(paste0("# interactive: ", interactive(), "\n"))
cat(paste0("# commandline: ", paste0(commandArgs(TRUE), collapse=" "), "\n"))
cat(paste0("# date: ", Sys.time(), "\n"))

# output the header.
cat(sprintf("gene\traw\trpm\tarea\n"))
for (genename in genes) {
    mat <- depthFromBam(genename, gtf, files, aVec)
    if (is.null(mat))
        next
    if (any(is.na(mat))) stop("read depth matrix contains NAs!   Stopping.")
    nreads <- attr(mat, "nreads")

    # try various normalization methods.
    # http://evolutionontheway.wordpress.com/2013/07/24/a-comprehensive-evaluation-of-normalization-methods-for-illumina-high-throughput-rna-sequencing-data-analysis/
    #
    pvalue.raw <- qarp.pvalue(mat, aVec)
    
    mat.rpm <- scale(mat, center=FALSE, scale=totals/1e6)
    pvalue.rpm <- qarp.pvalue(mat.rpm, aVec)

    # scale by the number of mapped reads in the gene
    area <- apply(mat, 2, function(r) sum(r))
    area[area==0] <- 1
    mat.area <- scale(mat, center=FALSE, scale=area)
    pvalue.area <- qarp.pvalue(mat.area, aVec)

    if (any(nreads == 0))
      genename <- paste0(genename, "*")
    
    cat(sprintf("%s\t%0.4g\t%0.4g\t%0.4g\n", genename, pvalue.raw, pvalue.rpm, pvalue.area))

    ## qarp.plotProfile(mat.rpm, aVec, showsplices=TRUE, invert=TRUE, 
    ##                  main=sprintf("SOC chemoresistive vs. chemosensitive\npvalue = %f", pvalue.rpm),
    ##                  xlab=paste0("position on ", genename),
    ##                  ylab="depth of reads (RPM)")
    ## qarp.plotMDS(mat.rpm, aVec)
    ## # annotate the graph with the pvalue for this gene
    ## mtext(sprintf("pvalue = %0.4g", pvalue.rpm), side=1, line=1, adj=0)
    ## mtext(sprintf("%s", genename), side=1, line=1, adj=1)

}
