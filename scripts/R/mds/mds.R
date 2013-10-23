require(vegan)
require(reshape)   # for melt
require(ggplot2)
require(grid)  # for "cm" units used when plotting



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
qarp.pvalue <- function(mat, aVec, method=NULL, permutations=NULL) {
    
    distance <- qarp.distance(mat, method=method)   
 
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

    if (is.null(permutations))
        df <-  adonis(distance ~ aVec)
    else {
        df <-  adonis(distance ~ aVec, permutations=permutations)
    }
    pvalue <- df$aov.tab[6][1,1]

    return(pvalue)
}
##' Calculate a distance matrix for a collection of points.
##'
##' At present, this routine is a thin veneer on the 'dist' function.  The indirection is useful only because we may want to replace the distance function at
##' some point so it will be conventient if all distance calculations go through one function.
##' @param mat matrix of points, each row is a single point and each column is a single dimensional component of that point.
##' @param method 
##' @return 
##' @author Chris Warth
qarp.distance <- function(mat, method=NULL) {
    if (is.null(method))
        distance <- dist(mat, method="euclid")
    else
        distance <- dist(mat)
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
qarp.read.cachedtsv <- function(filename, sep='\t', cache=TRUE) {
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
        ## HACK - remove the two rows of totals from the dataset....
        data <- data[c(-1,-2),]
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


##' Calculate pvalues for collections of RNA/RPT profiles.
##'
##' if genelist is not null, calculates pvalues only for the genes in
##' the list.  Otherwise pvalues will be calculated for all genes for
##' which there arwe entries in data.
##' 
##' The data frame containing per-gene read depth data.  Each row in
##' the data frame should have a \code{symbol} column and a
##' \code{tranPos} column indicating a gene and transcript position.
##' The combination of \code{symbol} and \code{tranPos} should be
##' unique in the data file.  Other columns give the number of reads
##' overlapping position \code{transPos} in the gene transcript.
##' Missing data rows are assumed have no reads at that position.
##'
##' @param data 
##' @param genelist
##' @return 
##' @author Chris Warth
qarp <- function(data, aVec, genelist=NULL, permutations=NULL, meancenter=TRUE) {
    if (is.null(genelist)) {
        ## extract the list of genes.
        genelist = levels(data$symbol)
        stopifnot(!is.null(genelist))
    }
    
    pvalues = sapply(X = genelist,
        FUN = function(gene, data, aVec, permutations) {
            mat <- qarp.matrix(data, gene, aVec, meancenter=meancenter)
            qarp.pvalue(mat, aVec, permutations=permutations)
        }, data = data, aVec = aVec, permutations = permutations)
    data.frame(row.names=genelist, pvalues)
}

##' Create a matrix suitable 
##' from a larger dataset.   The resulting matrix is suitable for
##' use in other QARP functions.
##'
##' The resulting matrix will be restricted to data for a single gene.  
##' 
##' @param data a data frame containing raw read counts for multiple positions along many genes.
##' @param gene the common name of a gene found in the data.
##' @param aVec a list of factors corresponsing to datasets in the
##' data frame.  This is used to identify wwhich columns in data
##' correspond to treated and control datasets.  There should be only
##' two factor levels in the list and each entry should have a name
##' corresponding to one of the name of one column in the data frame.
##' @param meancenter TRUE if the read count data for each gene should should be centered around the mean of the counts for that gene.
##' @return Returns a matrix suitable for use in the other QARP functions like qarp.plotProfile.
##' @author Chris Warth
qarp.matrix <- function(data, gene, aVec, meancenter=TRUE) {
    stopifnot(all(c("symbol", "transPos", "exonNumber") %in% names(data)))

    data.gene <- data[data$symbol==gene,]
    data.gene <- data.gene[order(data.gene$transPos),]		# order by transcript position
    rownames(data.gene) = data.gene$transPos
    isplices <- which(diff(data.gene$exonNumber) != 0)
    splices <- data[isplices, 'transPos']
    data.gene = data.gene[, names(aVec)]
    mat = as.matrix(t(data.gene))
    ## center read depth about the mean by default.
    if (meancenter)
        mat <- mat + 1 # add a pseudo count
        mat <- t(scale(t(mat), center=TRUE, scale=TRUE))
    mat[is.na(mat)] <- 0
    
    ## attach attributes to the matrix
    attr(mat, "gene") <- gene
    attr(mat, "filename") <- attr(data, "filename")
    attr(mat, "junctions") <- splices
    
    return(mat)
}


plotPvalue <- function(filename, datasets) {
    ## read the dataset from a file
    message("reading dataset...", appendLF = FALSE)
    data <- qarp.read.cachedtsv(filename)
    message("done.")

    # check that the dataset has some required columns
    stopifnot(all(c("symbol", "transPos", "exonNumber") %in% names(data)))

    # eliminate extra columns that we don't need.
    data = data[, c("symbol", "transPos", "exonNumber", datasets)]
    data$symbol <- as.factor(data$symbol)
    
    pvalues <- qarp(data, aVec)
 
    hist(pvalues, n=50, main=paste(filename))
    Sys.sleep(0.001)
    pvalues <- pvalues[order(pvalues)]
    data.frame(row.names=NULL, gene=names(pvalues), pvalues)
}


##' plot the read profiles across a gene.
##'
##' .. content for \details{} ..
##' @param mat read depths for a single gene.
##' @param aVec 
##' @param xlim limitations of x-axis, if applicable.   This allows you to zoom in on particular areas of interest.
##' @return Does not return anything of value.
##' @author Chris Warth
qarp.plotProfile <- function(mat, aVec, showsplices=TRUE, xlim=NULL) {
    gene <- attr(mat, "gene")
    colorPalette <- c("red", "blue")

    gg <- ggplot(melt(mat), aes(name="", x=X2,
                                y=value,
                                group=X1), environment = environment())
    gg <- gg + theme_bw()
    gg <- gg + theme(legend.title = element_text(size = 16, face = "bold"),
                     legend.text = element_text(size = 14, face = "bold"),
                     legend.position="top",
                     legend.direction="horizontal")
    gg <- gg + theme(legend.key = element_rect(size = 0.5, linetype="blank"))
    gg <- gg + theme(panel.border = element_blank())
    gg <- gg + theme(plot.margin = unit(c(0,0,0,0), "cm"))

    gg <- gg + ylab("read depth (normalized to RPM)")
    gg <- gg + xlab("Transcript position")
    gg <- gg + ylim(-max(mat),max(mat))
    gg <- gg + scale_colour_manual(name=gene, values=colorPalette,
                                    labels=c("Control", "Treated"))
    gg <- gg + geom_line(data=melt(mat[treated,,drop=FALSE]),aes(color=colorPalette[[1]]))
    gg <- gg + geom_line(data=melt(-mat[control,,drop=FALSE]), aes(color=colorPalette[[2]]))
    if (!is.null(xlim))
        gg <- gg + coord_cartesian(xlim=xlim)

    if (showsplices) {
        splices <- attr(mat, "junctions")
        if (!is.null(splices) && length(splices) > 0) {
            gg <- gg + geom_vline(xintercept = splices, alpha=.25,  color="red", linetype="dashed")
        }
    }
    print(gg)


}


##' plot multi-dimensional sampling plot.
##'
##' .. content for \details{} ..
##' @param mat 
##' @param aVec 
##' @return 
##' @author Chris Warth
qarp.plotMDS <- function(mat, aVec) {
    distance <- qarp.distance(mat)
    
    # get the two dismension points to be drawn on the scatter plot
    fit <- cmdscale(distance, eig=TRUE, k=2)
    df <- data.frame(fit$points)

    palette <- c("red", "blue")

    # make a factor colum that indicates the condition (treated or
    # control) of each dataset
    cond <- aVec

    # add a column indicating what color should be used for each point.
    df <- transform(df, col=palette[cond])

    par(mar=c(3, 3, 0.5, 1))  # Trim margin around plot [bottom, left, top, right]

    par(mgp=c(1.5, 0.2, 0))  # Set margin lines; default c(3, 1, 0) [title,labels,line]
    par(xaxs="r", yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

    plot(X2 ~ X1, data=df, col=as.character(col), pch=21, bg=as.character(col),
         xlab="", ylab="", cex=1.5, frame.plot=F, yaxt="n")

    ticks = pretty(c(min(df$X2), max(df$X2)), 3)

    # par(mgp=c(axis.title.position, axis.label.position, axis.line.position))
    mgp <- par("mgp")
    mgp[2] <- 0.5
    par(mgp=mgp)

    axis(2, at=ticks, labels=T, lwd=0, lwd.ticks=1, lty="solid", las=1, cex.axis=0.9)
    grid()

    # annotate the graph with the pvalue for this gene
    pvalue <- qarp.pvalue(mat, aVec)
    mtext(sprintf("pvalue = %0.4g", pvalue), side=1, line=1, adj=0)
}

