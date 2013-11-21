# load our libraries 
suppressMessages(library(IRanges))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))
source("~/Morris-Lab/scripts/R/mds/R/mds.R")


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

##' Plot a single sollection of depth measures across a gene.   This is not the routine to use to plot
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
calcBamDepth <- function(bamfile, gene, drop=NA) {
    stopifnot(is(gene, "GRanges"))
    stopifnot(is(bamfile, "character"))

    exons <- gene
    if (!is.na(drop) & drop <= length(exons))
      exons <- exons[-c(drop)]

    ## NOTE - RNA-seq alignments are often stranded, but the strand is
    ## not biologically significant.  Strand of a given RNA-seq
    ## fragment is a function of which round of PCR amplification was
    ## completed before sequencing.
    ##
    ## To select only those alignments that are in the same direction
    ## as the gene annotation, add the following to the scanBamFlag
    ## parameter of ScanBamParam :
    ##
    ##   isMinusStrand = (unique(as.vector(strand(gene))) == "-")
    
    param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
                          which = exons,
                          flag = scanBamFlag(isUnmappedQuery = FALSE))
    reads <- readGAlignmentsFromBam(bamfile, param=param)
    depth <- calcDepth(reads, gene)

    attr(depth, "nreads") <- length(reads)
    return(depth)
}

depthFromBam <- function(genename, gtf, files) {
    gene <- gtf[mcols(gtf)$gene_id==genename,]
    if (length(gene) == 0) {
        warning("gene \"", genename, "\" not found in list of features table 'gtf'.")
        return()
    }
    depths <- lapply(files, calcBamDepth, gene)
    # Here sapply() will return a matrix with samples across the columns,
    # and positions down the rows.
    mat <- sapply(depths, function(d) Reduce(c, sapply(d, as.vector)))
    attr(mat, "nreads") <- sapply(depths, function(d) attr(d, "nreads"))
    attr(mat, "junctions") <- junctions(gene)
    attr(mat, "gene") <- genename
    return(mat)
}
