# load our libraries 
library(IRanges) 
library(GenomicRanges) 
library(Rsamtools)
library(rtracklayer)
source("~/Morris-Lab/scripts/R/mds/R/mds.R")
## library(qarp) 

options(max.print=30)

## Holy !*#*@*, Rsamtools::bamCount is slow!   I actually don't know how long it takes to count the
## number of alignments in a decently sized BAM file because I got tired of waiting.
## Here is a function for doing the same thing in a few seconds.  
fastBamCount <- function(file) {
    cmd<-paste0("samtools idxstats ", file, "| awk '{ sum += $3+$4 } END {print sum}'")
    cnt <- as.numeric(system(cmd, intern=TRUE))
}

junctions <- function(gene) {
    stopifnot(is(gene, "GRanges"))
    j <- Reduce(sum, width(gene), accumulate=T)
    length(j) <- length(j)-1
    # genomic numbering is 1-based, so shift all junctions right by one
    return(j+1)
}

## insert a new row into a data frame at a particular row location.
## See http://stackoverflow.com/a/11587051/1135316
##
insertRow <- function(existingDF, newrow, r) {
    existingDF <- rbind(existingDF,newrow)
    existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
    row.names(existingDF) <- 1:nrow(existingDF)
    return(existingDF)  
}

## A paradigm for collapsing the list-of-lists into a single list is

.unlist <- function (x)
{
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)) {
        structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
        do.call(c, x)
    }
}
##' Plot a single sollection of depth mesures across a gene.   This is not the routine to use to plot
##' Multiple datasets of depth across a gene...
##'
##' .. content for \details{} ..
##' @title 
##' @param depth 
##' @param showJunctions 
##' @return 
##' @author Christopher Warth
plotDepth <- function(depth, showJunctions=TRUE) {
    vdepth <- Reduce(c, sapply(depth, as.vector))
    plot(vdepth, type="l")
    if (showJunctions) {
        j <- Reduce(sum, sapply(depth, function(d) sum(runLength(d))), accumulate=T)
        length(j) <- length(j)-1
        j <- j+1
        abline(v=j, col="red")
    }
}

calcDepth <- function(algns, gene) {
    stopifnot(is(algns, "GAlignments"))
    stopifnot(is(gene, "GRanges"))

    exons <- ranges(gene)
#    readranges <- sapply(reads, function(x) IRanges(x$pos, width=x$qwidth))
    readranges <-  IRanges(start(algns), width=width(algns))
    # calculate the coverage for one contiguous range.
    # the coverage is shifted so it begins at 0.
    calcCoverage <- function(iexon, rr) {
         coverage(rr, shift=-start(exons[iexon])+1, width=width(exons[iexon]))
    }

   depth <- sapply(seq_along(exons), calcCoverage, rr=readranges)
}

##' Calculate the depth of coverage across all a set of features.  USually this means calculating
##' dpeth of coverage over a set of exons, but other features that can be represented in GRanges objects can also be used.
##'
##' Returns a list of RLE (Run-Length Encoded) structures, one per exon, representing the depth of coverage
##' across each of the features, usually across exons of a gene.
##' @title 
##' @param bamfile character name of a bam file, e.g. "~/data/foo.bam".   There should be a corresponding index file in the same location, e.g. "~/data/foo.bam/bai"
##' @param gene a GRanges object specifying the genomic ranges of features over which coverage depth should be calculated.
##' Usually these are read from a GTF file with \code{rtracklayer::import()} and a subset of ranges specific to a particular gene are passed in to this routine.
##' @return Returns a list of RLE (Run-Length Encoded) structures, one per subrange,
##' representing the depth of coverage across each of the exons of the
##' gene.  See
##' \code{\url{http://www.bioconductor.org/packages/release/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf}}
##' for more information about the RLE structure
##' @author Christopher Warth
calcBamDepth <- function(bamfile, gene) {
    stopifnot(is(gene, "GRanges"))
    stopifnot(is(bamfile, "character"))
    param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
                          which = gene,
                          flag = scanBamFlag(isUnmappedQuery = FALSE,
                                             isMinusStrand = (unique(as.vector(strand(gene))) == "-")
                            ))
#    reads <- scanBam(bamfile, param = param)
    reads <- readGAlignmentsFromBam(bamfile, param=param)
    depth <- calcDepth(reads, gene)
}


## read SOC clinical outcomes spreadsheet
prefix <- file.path("", "Volumes", "homes")
print(sprintf("prefix = %s", prefix))
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


genename <- "ACTB"
gtf <- rtracklayer::import(gtfname, genome="hg19", feature.type=c("exon"), asRangedData=FALSE, format="gtf")
# delete non-canonical chromosomes (haplotypes and such)
# only keep features whose chromosome labels don't have underscores.
# e.g.  we throw out "chr3_hap6" but keep "chr3"
# Note: we should probably get rid of chrM, mitochondrial features, as well.
gtf <- gtf[seqnames(gtf) %in% levels(seqnames(gtf))[grep("_|chrM", levels(seqnames(gtf)), invert=TRUE)]]

pvalueFromBam <- function(genename, files, aVec) {
    gene <- gtf[mcols(gtf)$gene_id==genename,]
    depths <- lapply(files, calcBamDepth, gene)
    totals <- sapply(files, fastBamCount)
    
    # Here sapply() will return a matrix with samples across the columns,
    # and positions down the rows.
    mat <- sapply(depths, function(d) Reduce(c, sapply(d, as.vector)))
    mat <- scale(mat, center=FALSE, scale=totals)
    

    attr(mat, "junctions") <- junctions(gene)
    attr(mat, "gene") <- genename

    ## qarp.plotProfile(mat, aVec, showsplices=TRUE, invert=TRUE, 
    ##                  main="SOC chemoresistive vs. chemosensitive", xlab=paste0("position on ", genename),
    ##                  ylab="depth of reads\nnormalized to total alignments (RPM)")
    
    qarp.pvalue(mat, aVec)
}

pvalueFromBam("ACTB", files, aVec)
pvalueFromBam("Pbsn", files, aVec)


