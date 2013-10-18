require(vegan)
require(xlsx)



##' Calculate a psuedo p-value for a set of categorized RNA-SEQ profiles.
##'
##' This routine calculates a single measure indicating the relative
##' similarity of RNA-SEQ profiles.  Profiles are categorized
##' according to factors in \Code{aVec}.  If the profiles within each
##' category are more similar to others in the same category than to
##' profiles in other categories, then the calculated value will be
##' relatively high.  If the profiles in one group are self-similar
##' but distinct from other groups then the calculated pvalue will be
##' low.  This suggests the profiles are affected but the experimental
##' conditions.
##' 
##' @param mat
##' 
##' numeric matrix of read depth at each position along a transcript.
##' Each row is an experiment, and each column is a transcript position.  
##' 
##' @param aVec
##' 
##' vector of factors that categorize the experiments.  The order of
##' elements in \Code{aVec} should be the same as the order of rows in
##' \Code{mat}.  Typically there will be two factors that divide the
##' available experiments into two groups, e.g. \Code{"treated",
##' "control"} .  More factors can be used but nested factors are not
##' supported.
##' 
##' @param method
##'
##' the distance measure to be used. This must be one of "euclidean",
##' "maximum", "manhattan", "canberra", "binary" or "minkowski". Any
##' unambiguous substring can be given.
##' 
##' @param ...
##' remaining parameters are passed to vegan::adonis
##' @return 
##' @author Chris Warth
##' @export
qart.pvalue <- function(mat, aVec, method="euclid", permutations=10000, ...) {
    
    distance <- dist(mat, method=method)
    aVec <- as.factor(sapply(rownames(mat) %in% control, function(x) if (x) {"control"} else {"treated"}))

    # There are two ways to call vegan::adonis.  The first used raw data
    # and the second uses a distance matrix calculated from the raw data.
    # In the former case, you pass a matrix of raw measurments to
    # adonis() and it calculates the distance matrix for you using the
    # algorithm specified by the 'method' argument.  In the second type
    # of usage, you precompute a distance matrix from the raw data using
    # any algorithm you like and then pass the distance matrix as the
    # first argument to adonis().  adonis() can tell which type of usage
    # you expect by looking at the type of the first argument; if it is
    # of class 'dist' then you have precomputed a distance matrix,
    # otherwise you want adonis() to do it for you.
    #
    # Why would you want to use adonis() one way versus the other?  It is
    # usually going to be easier to pass raw data to adonis() and let it
    # calculate the distance matrix for you.  However, if you want to try
    # out a distance calculation that is not one of the algorithms
    # available through adonis() then you will want to precopute your own
    # distance matrix and pass it in.

    df <-  adonis(distance ~ aVec, permutations=permutations, ...)
    pvalue <- df$aov.tab[6][1,1]

    return(pvalue)
}

##'  Read a tab-seperated-value file and cache the resulting dataset so reading it will be faster next time.
##'
##' .. content for \details{} ..
##' @param the name of the file which the data are to be read from.   See \Code{help('read.csv')} for more information.
##' @param sep the field separator character.  Default is a single tab character.  See \Code{help('read.csv')} for more information.
##' @param cache if TRUE, attempt to read and write a cached version of the file in the same
##' directory; otherwise no cache is used.
##' @return
##' dataframe with the data from the data file.
##' @author Chris Warth
##' @export
qart.read.cachedtsv <- function(filename, sep='\t', cache=TRUE) {
    data <- NULL
    base <- sub("^(.*)\\..*", "\\1", filename)
    rdafile = paste0(base, ".Rda")
    
    # try to read the ached version
    if (cache) {
        e <- try({load(paste0(base, ".Rda"))}, silent=TRUE)
        if (class(e) != "try-error") {
            attr(data, "filename") = filename
            attr(data, "cachefile") = rdafile
            attr(data, "readcache") = TRUE
        }
    }
    if (is.null(data)) {
        data <- read.csv(filename, sep=sep, stringsAsFactors=FALSE)
        attr(data, "filename") = filename
        if (cache) {
            attr(data, "cachefile") = rdafile
            # if we can, save a copy in .rda format.
            if (file.access(rdafile, 2)) {
                save("data", file=rdafile)
                attr(data, "wrotecache") = TRUE
            } else {
                warning("Cannot cache .Rda file")
                attr(data, "wrotecache") = FALSE
            }
        }
    }

    return(data)
}

##' Calculate a psuedo p-value for a single gene within a larger dataset.
##'
##' Use a pseudo f-test to estimate a p-value statistic categorized
##' points in the distance matrix.  A null distribution is constrcted
##' by randomly permuting the categories associated with each data
##' point in the distance matrix.  By mixing the categories and data
##' points many times, a null distribution of distance matricies
##' can be estimated.
##' 
##' @param gene
##' 
##' Identifier in \code{data$symbol} for read depth data for a single
##' gene.  Usually this is a character string, less often a factor.
##' 
##' @param data
##' 
##' The data frame containing per-gene read depth data.  Each row in
##' the data frame should have a \code{symbol} column and a
##' \code{tranPos} column indicating a gene and transcript position.
##' The combination of \code{symbol} and \code{tranPos} should be
##' unique in the data file.  Other columns give the number of reads
##' overlapping position \code{transPos} in the gene transcript.
##' Missing data rows are assumed have no reads at that position.
##' 
##' @return
##' 
##' Returns a pseudo pvalue estimating the probability of seeing a
##' similar categorization pattern of read profiles by random chance.
##' 
##' @author Chris Warth
calcPvalue <- function(gene, data = NULL) {
    ## Make a data matrix, where each row is an experiment, and each
    ## column is a transcript position.  The value at each position
    ## represents the read depth at that particular position (col) in one
    ## experiment (row).
    data.gene <- data[data$symbol==gene,]
    data.gene <- data.gene[order(data.gene$transPos),]		# order by transcript position
    rownames(data.gene) = data.gene$transPos
    data = subset(data.gene, select=c(-symbol, -transPos))
    mat = as.matrix(t(data))


    aVec <- as.factor(sapply(rownames(mat) %in% control, function(x) if (x) {"control"} else {"treated"}))
    p <- qart.pvalue(mat, aVec)
    return(p)
}


plotPvalue <- function(filename, datasets) {
    ## read the dataset from a file
    message("reading dataset...", appendLF = FALSE)
    data <- qart.read.cachedtsv(filename)
    message("done.")

    # check that the dataset has some required columns
    stopifnot(all(c("symbol", "transPos") %in% names(data)))

    data = data[, c("symbol", "transPos", datasets)]
    data$symbol = as.factor(data$symbol)
    
    ## extract the list of genes.
    genes = levels(data$symbol)
    stopifnot(!is.null(genes))

    pvalues = sapply(X = genes, FUN = calcPvalue, data = data)
 
    hist(pvalues, n=50, main=paste(filename))
    Sys.sleep(0.001)
    pvalues <- pvalues[order(pvalues)]
    data.frame(row.names=NULL, gene=names(pvalues), pvalues)
}



filenames = c("PEO1-RPT-corrUp-coverage.tsv", "PEO1-RPT-corrDown-coverage.tsv", "PEO1-RPT-top200-coverage.tsv")

## specify which datasets are grouped together
control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
datasets = c(control = control, treated = treated)



pl <- lapply(filenames, plotPvalue, datasets)
pvalues.best <- do.call(cbind, as.list(lapply(pl, function(df) head(df,n=10))))
pvalues.worst <- do.call(cbind, as.list(lapply(pl, function(df) { df <- tail(df, n=10)
                                                                  df <- df[order(df$pvalue, decreasing=TRUE),]})))

write.xlsx(pvalues.best, file = "pvalues.best.xlsx", sheetName = "best", row.names = FALSE)
write.xlsx(pvalues.worst, file = "pvalues.worst.xlsx", sheetName = "worst", row.names = FALSE)


