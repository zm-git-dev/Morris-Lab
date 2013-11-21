##' Calculates pvalues for RNA-SEQ profiles from chemosensitive and chemoresistive phenotypes of ovarian tumor.
##'
##' @usage   R CMD BATCH --slave chemotype_ytest.R output.txt
##' @author  Chris Warth
##' @date 11-11-2013


source("qarputils.R")

## read SOC clinical outcomes spreadsheet
prefix <- file.path("", "Volumes", "homes")
if (!file.exists(file.path(prefix, "data"))) {
    prefix <- file.path("~")
}
file.survival = file.path(prefix, "data", "SOC", "2013-09-25_McIntoshSOC_Survival.csv")
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

genes <- c("ACTB", "GAPDH", "MALAT1", "EEF1A1", "RPLP0", "RPL41", "EEF1G",
          "CD74", "AHNAK")
genes <- c("CNTNAP2:MIR548F3:MIR548T", "QTRTD1", "WBSCR16")

print(genes)

totals <- sapply(files, fastBamCount)

# set the random seed to experiments are reproducible
set.seed(1234)

# output the header.
cat(paste0("# BAM files: ", paste0(files, collapse=", "), "\n"))
cat(paste0("# phenotypes: ", paste0(aVec, collapse=", "), "\n"))
cat(paste0("# annotations: ", gtfname, "\n"))
cat(paste0("# gene count: ", length(genes), "\n"))
cat(paste0("# total alignments: ", paste0(totals, collapse=", "), "\n"))
cat(paste0("# date: ", Sys.time(), "\n"))
# cat(paste0("# script: ", sys.frame(1)$ofile, "\n"))

# output the header.
cat(sprintf("gene\traw\trpm\tarea\n"))
for (genename in genes) {
    mat <- depthFromBam(genename, gtf, files, aVec)
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

    par(mfrow=c(2,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .1, 2.5, .1), cex=1, las=1)
    qarp.plotProfile(mat, aVec, showsplices=TRUE, invert=TRUE,
                     main="", xlab="", ylab="", yaxt="n", xaxt="n"
                     ## main=sprintf("SOC chemoresistive vs. chemosensitive\npvalue = %f", pvalue),
                     ## xlab=paste0("position on ", genename),
                     ## ylab="depth of reads (RAW)"
                     )
    mtext(sprintf("gene: %s", genename), side=3, line=0, adj=0)
    mtext(sprintf("normalized area"), side=2, line=0, rot=90)

    qarp.plotMDS(mat, aVec,  main="", xlab="", ylab="", yaxt="n", xaxt="n")
    mtext(sprintf("pvalue: %0.4g", pvalue.raw), side=3, line=0, adj=1)

    qarp.plotProfile(mat.area, aVec, showsplices=TRUE, invert=TRUE,
                     main="", xlab="", ylab="", yaxt="n", xaxt="n"
                     ## main=sprintf("SOC chemoresistive vs. chemosensitive\npvalue = %f", pvalue),
                     ## xlab=paste0("position on ", genename),
                     ## ylab="depth of reads (RAW)"
                     )
    mtext(sprintf("gene: %s", genename), side=3, line=0, adj=0)

    qarp.plotMDS(mat.area, aVec,  main="", xlab="", ylab="", yaxt="n", xaxt="n")
    mtext(sprintf("pvalue: %0.4g", pvalue.area), side=3, line=0, adj=1)

}



