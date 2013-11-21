##' Draw a graphs coparing raw read profiles with profiles normalized by area.
##' 
##' @usage   R CMD BATCH --slave --quiet mkpics.R
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

totals <- sapply(files, fastBamCount)

# now read the GTF file to get descriptions of gene regions
gtfname <- system.file("python", "hg19-refSeq-union.gtf", package="qarp")
stopifnot(nchar(gtfname)>0)
gtf <- rtracklayer::import(gtfname, genome="hg19", feature.type=c("exon"), asRangedData=FALSE, format="gtf")
# delete non-canonical chromosomes (haplotypes and such)
# only keep features whose chromosome labels don't have underscores.
# e.g.  we throw out "chr3_hap6" but keep "chr3"
# Note: we should probably get rid of chrM, mitochondrial features, as well.
gtf <- gtf[seqnames(gtf) %in% levels(seqnames(gtf))[grep("_|chrM", levels(seqnames(gtf)), invert=TRUE)]]

samples <- read.csv("top300_1.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, comment.char="#")
samples[,2] = as.numeric(samples[,2])
samples[,3] = as.numeric(samples[,3])
samples[,4] = as.numeric(samples[,4])

genename <- "WFDC2"
png(sprintf("%s.png", genename), width = 800, height = 480)

set.seed(1234)   # make the experiments repeatable, as long as permutations are the same.
mat <- depthFromBam(genename, gtf, files, aVec)
if (any(is.na(mat))) stop("read depth matrix contains NAs!   Stopping.")
nreads <- attr(mat, "nreads")

par(mfrow=c(1,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .1, 2.5, .1), cex=1, las=1)
## draw grahs based on raw read count
qarp.plotProfile(mat, aVec, showsplices=TRUE, invert=TRUE,
                 main="", xlab="", ylab="", yaxt="n", xaxt="n"
                 ## main=sprintf("SOC chemoresistive vs. chemosensitive\npvalue = %f", pvalue),
                 ## xlab=paste0("position on ", genename),
                 ## ylab="depth of reads (RAW)"
                 )
mtext(sprintf("gene: %s", genename), side=2, line=0, las=3, cex=1.5)
mtext(sprintf("raw"), side=3, line=0)

# scale by the number of mapped reads in the gene
area <- apply(mat, 2, function(r) sum(r))
area[area==0] <- 1
mat.area <- scale(mat, center=FALSE, scale=area)
pvalue.area <- qarp.pvalue(mat.area, aVec)

## draw grahs based on raw read count
qarp.plotProfile(mat.area, aVec, showsplices=TRUE, invert=TRUE,
                 main="", xlab="", ylab="", yaxt="n", xaxt="n"
                 ## main=sprintf("SOC chemoresistive vs. chemosensitive\npvalue = %f", pvalue),
                 ## xlab=paste0("position on ", genename),
                 ## ylab="depth of reads (RAW)"
                 )
mtext(sprintf("area normalized"), side=3, line=0)

dev.off()
