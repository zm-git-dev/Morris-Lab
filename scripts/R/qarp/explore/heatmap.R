##' Produce a heatmap of intermediate result distance matrix used to compute pvalues.
##' 
##'
##' @usage   R CMD BATCH --slave --quiet heatmap.R
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

        ## draw grahs based on raw read count
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

        d.raw <- dist(t(mat))
        heatmap.2(as.matrix(d.raw), Rowv=FALSE, Colv=FALSE,
                  main=sprintf("%s raw data", genename),
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

        ## NOTE: heatmap and heatmap.2 call plot.new() internally, invalidating any layout directives you may have
        ## assembled prior to this point.   It is impossible with show two heatmaps next to one another.
        ## http://stackoverflow.com/questions/13081310/combining-multiple-complex-plots-as-panels-in-a-single-figure
        d.area <- dist(t(mat.area))
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
samples[,2] = as.numeric(samples[,2])
samples[,3] = as.numeric(samples[,3])
samples[,4] = as.numeric(samples[,4])

## draw the histograms

pdf("heatmap.pdf")

# set the graph margins so we can fit four graphs on a page,
# two-by-two, and have room for a title at the top of each pair of
# graphs.
#
# mar = c(bottom, left, top, right) 
par(mfrow=c(2,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .4, 2.5, .1), cex=1, las=1)

makeGraphs(c("WFDC2"))
dev.off()
