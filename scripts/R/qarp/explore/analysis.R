##' Analyze a list of pre-computed pvalues calculated from RNA/RIBO-SEQ profiles.
##' The idea is to first sort the pvalues of the profiles to find which genes look interesting and then
##' plot some of the results.
##'
##' @usage   R CMD BATCH --slave --quiet analysis.R
##' @author  Chris Warth
##' @date 11-11-2013
source("qarputils.R")

# Draw profile graphs and scatter plots for the genes in genelist
makeGraphs <- function(genes) {
    # set the graph margins so we can fit four graphs on a page,
    # two-by-two, and have room for a title at the top of each pair of
    # graphs.
    #
    # mar = c(bottom, left, top, right) 
    par(mfrow=c(2,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .1, 2.5, .1), cex=1, las=1)
    for (genename in genes) {
        set.seed(1234)   # make the experiments repeatable, as long as permutations are the same.
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

        ## draw graphs based on raw read count
        qarp.plotProfile(mat, aVec, showsplices=TRUE, invert=TRUE,
                         main="", xlab="", ylab="", yaxt="n", xaxt="n"
                         ## main=sprintf("SOC chemoresistive vs. chemosensitive\npvalue = %f", pvalue),
                         ## xlab=paste0("position on ", genename),
                         ## ylab="depth of reads (RAW)"
                         )
        mtext(sprintf("gene: %s", genename), side=3, line=0, adj=0)
        mtext(sprintf("RAW"), side=2, line=0, las=3, cex=1.5)

        qarp.plotMDS(mat, aVec,  main="", xlab="", ylab="", yaxt="n", xaxt="n")
        mtext(sprintf("pvalue: %0.4g", pvalue.raw), side=3, line=0, adj=1)

        ## draw graphs based on read counts normalized to the area under the graph
        qarp.plotProfile(mat.area, aVec, showsplices=TRUE, invert=TRUE,
                         main="", xlab="", ylab="", yaxt="n", xaxt="n"
                         ## main=sprintf("SOC chemoresistive vs. chemosensitive\npvalue = %f", pvalue),
                         ## xlab=paste0("position on ", genename),
                         ## ylab="depth of reads (RAW)"
                         )
        mtext(sprintf("gene: %s", genename), side=3, line=0, adj=0)
        mtext(sprintf("normalized to area"), side=2, line=0, las=3, cex=1.5)

        qarp.plotMDS(mat.area, aVec,  main="", xlab="", ylab="", yaxt="n", xaxt="n")
        mtext(sprintf("pvalue: %0.4g", pvalue.area), side=3, line=0, adj=1)
    }
    cat("\n")
}


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
# some genes will have '*' in the name if some experiments have no reads for that gene
# strip and ignore this flag for now.
samples[,1] <- sub("*", "", samples[,1], fixed=TRUE)
samples[,2] = as.numeric(samples[,2])
samples[,3] = as.numeric(samples[,3])
samples[,4] = as.numeric(samples[,4])

#pdf("analysis.pdf")
png(sprintf("hist.png"), width = 800, height = 480)

## draw the histograms

# set the graph margins so we can fit four graphs on a page,
# two-by-two, and have room for a title at the top of each pair of
# graphs.
#
# mar = c(bottom, left, top, right) 
par(mfrow=c(2,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .4, 2.5, .1), cex=1, las=1)

hist(samples[,2], main="")
abline(v=.05, col="red")
mtext(sprintf("raw read counts"), side=3, line=0)

hist(samples[,3], main="")
abline(v=.05, col="red")
mtext(sprintf("normalized read counts (%s)", names(samples)[3]), side=3, line=0)

hist(samples[,4], main="")
abline(v=.05, col="red")
mtext(sprintf("normalized read counts (%s)", names(samples)[4]), side=3, line=0)

# order the genes according to the pvalues derived from raw read counts.
samples <- samples[order(samples[,2]),]

# nsamples:  # of genes to sample from both the lowest and highest pvalues.
# if nsamples == 3, you would get 6 lines of output, 3 from high pvalues and 3 from low pvalues.
nsamples=3
genes <- unique(c(samples[1:nsamples,1], rev(samples[,1])[1:nsamples]))

makeGraphs(genes)

# order the genes according to the pvalues derived from read counts normalized by area.
samples <- samples[order(samples[,4]),]
genes <- unique(c(samples[1:nsamples,1], rev(samples[,1])[1:nsamples]))

makeGraphs(genes)


dev.off()


