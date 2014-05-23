#!/usr/bin/env Rscript
##' Show the effect of depletion on pvalues
##'
##' @usage   depletion [--help] [--genes <gene>,<gene2>,...] [--genefile <genefile>] [--sample <integer>]
##' @author  Chris Warth
##' @date 11-11-2013
library(proto)
suppressMessages(library(optparse))

# if executing locally, shared directories will be mounted under "/Volumes".
# if on the rhino servers, just go to the home directory.
prefix <- file.path("", "Volumes", "homes")
if (!file.exists(file.path(prefix, "data"))) {
    prefix <- file.path("~")
}

if (interactive()) {
    samplesize <- 5
    genes.file <- file.path(prefix, "data", "SOC", "genelist.txt")
    genes.file <- file.path(prefix, "data", "SOC", "isoforms.txt")
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
deplete <- function(genename, gtf, files, aVec) {
    mat <- depthFromBam(genename, gtf, files, aVec)
    if (is.null(mat))
        next
    if (any(is.na(mat))) stop("read depth matrix contains NAs!   Stopping.")

    start <- as.integer(nrow(mat)*.2)
    end <- as.integer(nrow(mat)*.4)
    datasets <- (aVec == levels(aVec)[1])
    for (depletion in c(0, .1, .2, .4, .6, .8)) {
        # remove reads from the middle of one of the datasets
        depleted.region <- mat[start:end,datasets]
        depleted.depth <- apply(depleted.region, 2, max)
        depleted.region <- as.integer(depleted.region * (1-depletion))
        
        mat.deplete <- mat
        mat.deplete[start:end, datasets] <- depleted.region
        area <- apply(mat.deplete, 2, function(r) sum(r))
        area[area==0] <- 1
        mat.area <- scale(mat.deplete, center=FALSE, scale=area)
        pvalue.area <- qarp.pvalue(mat.area, aVec)
        qarp.plotProfile(mat.area, aVec, showsplices=TRUE, invert=TRUE,
                         main="", xlab="", ylab="", yaxt="n", xaxt="n"
                         )
        if (depletion != 0) 
            mtext(sprintf("%d%% depleted", as.integer(depletion*100)), side=3, line=0, adj=.2, cex=.9)
        mtext(sprintf("pvalue %0.4g", pvalue.area), side=3, line=0, adj=.8, cex=.9)
    }
    mtext(sprintf("%s", genename), side=3, line=1, outer=TRUE)
    mtext(sprintf("pvalues with increasing levels of depletion", genename), side=3, line=0, outer=TRUE)

    mtext(sprintf("built with depletion.R"), side=1, line=1, outer=TRUE, col="lightgrey", cex=.7)
    
    
    ##cat(sprintf("%s\t%0.4g\t%0.4g\t%0.4g\n", genename, pvalue.raw, pvalue.rpm, pvalue.area))

    ## qarp.plotProfile(mat.rpm, aVec, showsplices=TRUE, invert=TRUE, 
    ##                  main=sprintf("SOC chemoresistive vs. chemosensitive\npvalue = %f", pvalue.rpm),
    ##                  xlab=paste0("position on ", genename),
    ##                  ylab="depth of reads (RPM)")
    ## qarp.plotMDS(mat.rpm, aVec)
    ## # annotate the graph with the pvalue for this gene
    ## mtext(sprintf("pvalue = %0.4g", pvalue.rpm), side=1, line=1, adj=0)
    ## mtext(sprintf("%s", genename), side=1, line=1, adj=1)

}

#png(sprintf("depletion.png"), width = 800, height = 1000)
par(mfrow=c(3,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .1, 2.5, .1), cex=1, las=1)
deplete("CNPY2:PAN2", gtf, files, aVec)
deplete("WFDC2", gtf, files, aVec)
#dev.off()

pdf("isoform.depletion.pdf")
par(mfrow=c(3,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .1, 2.5, .1), cex=1, las=1)
lapply(head(genes,n=10), deplete, gtf, files, aVec)
dev.off()
