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

draw <- function(genename, gtf, files, aVec) {
    mat <- depthFromBam(genename, gtf, files)
    if (is.null(mat))
        next
    if (any(is.na(mat))) stop("read depth matrix contains NAs!   Stopping.")

    area <- apply(mat, 2, function(r) sum(r))
    area[area==0] <- 1
    mat.area <- scale(mat, center=FALSE, scale=area)
    pvalue.area <- qarp.pvalue(mat.area, aVec)
    qarp.plotProfile(mat.area, aVec, showsplices=TRUE, invert=TRUE,
                     main="", xlab="", ylab="", yaxt="n", xaxt="n"
                     )
    mtext(sprintf("pvalue %0.4g", pvalue.area), side=3, line=0, adj=.8, cex=.9)
    mtext(sprintf("%s", genename), side=3, line=0, adj=.2, cex=.9)

    mtext(sprintf("built with slides.R"), side=1, line=1, outer=TRUE, col="lightgrey", cex=.7)
    
    
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

genes <- c("TLN1", "SF1", "IGFBP2", "WFDC2", "MUC16", "MSLN")

#pdf("isoform.depletion.pdf")
par(mfrow=c(3,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .1, 2.5, .1), cex=1, las=1)
lapply(head(genes,n=10), draw, gtf, files, aVec)
#dev.off()

png("WFDC2.png", height=3.79, width=8.58, units="in", res=78)
par(mar=c( .1, .1, 2.5, .1), cex=1, las=1)
genename <- "WFDC2"
mat <- depthFromBam(genename, gtf, files)
qarp.plotProfile(mat, aVec, showsplices=TRUE, invert=TRUE,
                 main="", xlab="", ylab="", yaxt="n", xaxt="n", bty="n"
                 )
#mtext(sprintf("pvalue %0.4g", pvalue.area), side=3, line=0, adj=.8, cex=.9)
#mtext(sprintf("%s", genename), side=3, line=0, adj=.2, cex=.9)
dev.off()


par(mfrow=c(2,2), oma=c(0, 0, 0, 0), mar=c( .1, .1, 2.5, 1), cex=1, las=1)
mat.rpm <- scale(mat, center=FALSE, scale=totals/1e6)
pvalue.rpm <- qarp.pvalue(mat.rpm, aVec)
qarp.plotProfile(mat.rpm, aVec, showsplices=TRUE, invert=FALSE,
                 main="", xlab="", ylab="", yaxt="n", xaxt="n", bty="n"
                 )
# mtext(sprintf("pvalue %0.4g", pvalue.rpm), side=3, line=0, adj=.8, cex=.9)
mtext(sprintf("%s", genename), side=3, line=0, adj=.2, cex=.9)
qarp.plotMDS(mat.rpm, aVec,  main="", xlab="", ylab="", yaxt="n", xaxt="n")


par(mfrow=c(2,2), oma=c(0, 0, 0, 0), mar=c( .1, .1, 2.5, 1), cex=1, las=1)
area <- apply(mat, 2, function(r) sum(r))
area[area==0] <- 1
mat.area <- scale(mat, center=FALSE, scale=area)
pvalue.area <- qarp.pvalue(mat.area, aVec)
qarp.plotProfile(mat.area, aVec, showsplices=TRUE, invert=TRUE,
                 main="", xlab="", ylab="", yaxt="n", xaxt="n", bty="n"
                 )
# mtext(sprintf("pvalue %0.4g", pvalue.rpm), side=3, line=0, adj=.8, cex=.9)
mtext(sprintf("%s", genename), side=3, line=0, adj=.2, cex=.9)
qarp.plotMDS(mat.area, aVec,  main="", xlab="", ylab="", yaxt="n", xaxt="n")

options(max.print=3)
d.raw <- dist(t(mat), method="euclid")
heatmap.2(as.matrix(d.raw), Rowv=FALSE, Colv=FALSE,
          main=sprintf("%s unnormalized", genename),
          dendrogram="none", trace="none", key=FALSE, scale="none",
          rowsep=c(10), colsep=c(10), 
          ## lmat=rbind(c(2),c(3),c(1),c(4)), 
          ## lhei=c(1,1,9,0), 
          ## lwid=c(1),
          add.expr=c(
              mtext("resistant", side=2, line=0, adj=.75, las=3, cex=1.5),
              mtext("sensitive", side=2, line=0, adj=.25, las=3, cex=1.5),
              mtext("sensitive", side=3, line=0, adj=.75, cex=1.5),
              mtext("resistant", side=3, line=0, adj=.25, cex=1.5)
              ))
mtext(genename, outer=TRUE)


d.area <- dist(t(mat.area), method="euclid")
heatmap.2(as.matrix(d.area), Rowv=FALSE, Colv=FALSE,
          main=sprintf("%s normalized by area", genename),
          dendrogram="none", trace="none", key=FALSE, scale="none",
          rowsep=c(10), colsep=c(10), 
          ## lmat=rbind(c(2),c(3),c(1),c(4)), 
          ## lhei=c(1,1,9,0), 
          ## lwid=c(1),
          add.expr=c(
              mtext("resistant", side=2, line=0, adj=.75, las=3, cex=1.5),
              mtext("sensitive", side=2, line=0, adj=.25, las=3, cex=1.5),
              mtext("sensitive", side=3, line=0, adj=.75, cex=1.5),
              mtext("resistant", side=3, line=0, adj=.25, cex=1.5)
              ))
mtext(genename, outer=TRUE)


png("cancer-specific-%03d.png", width=600, height=300)
par(mfrow=c(1,2), oma=c(0, 0, 0, 0), mar=c( .1, .1, 2.5, 1), cex=1, las=1)
for (genename in genes) {
    mat <- depthFromBam(genename, gtf, files)
    area <- apply(mat, 2, function(r) sum(r))
    area[area==0] <- 1
    mat.area <- scale(mat, center=FALSE, scale=area)
    pvalue.area <- qarp.pvalue(mat.area, aVec)
    qarp.plotProfile(mat.area, aVec, showsplices=TRUE, invert=TRUE,
                     main="", xlab="", ylab="", yaxt="n", xaxt="n"
                     )
    mtext(sprintf("pvalue %0.4g", pvalue.area), side=3, line=0, adj=.8, cex=1.2)
    mtext(sprintf("%s", genename), side=3, line=0, adj=.2, cex=1.2)
}
dev.off()
