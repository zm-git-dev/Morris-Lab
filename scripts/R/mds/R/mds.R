require(vegan)
require(reshape)   # for melt
require(grid)  # for "cm" units used when plotting
require(gplots)

## http://stackoverflow.com/a/12996160/1135316
addTrans <- function(color,trans)
  {
      # This function adds transparancy to a color.
      # Define transparancy with an integer between 0 and 255
      # 0 being fully transparant and 255 being fully visable
      # Works with either color and trans a vector of equal length,
      # or one of the two of length 1.

      if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
      if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
      if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

      num2hex <- function(x)
        {
            hex <- unlist(strsplit("0123456789ABCDEF",split=""))
            return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
        }
      rgb <- rbind(col2rgb(color),trans)
      res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
      return(res)
  }

##' Calculate a psuedo p-value for a set of categorized RNA-SEQ profiles.
##'
##' This routine calculates a single measure indicating the relative
##' similarity of RNA-SEQ profiles.  Profiles are categorized
##' according to factors in \code{aVec}.  If the profiles within each
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
##' elements in \code{aVec} should be the same as the order of rows in
##' \code{mat}.  Typically there will be two factors that divide the
##' available experiments into two groups, e.g. \code{"treated",
##' "control"} .  More factors can be used but nested factors are not
##' supported.
##' 
##' @param method
##'
##' the distance measure to be used. The documentation for the \code{distance} parameter in \code{qarp.distance}
##' 
##' @param ...
##' remaining parameters are passed to vegan::adonis
##' @return a floating point value representing the pvalue calculated for these profiles. 
##' @author Chris Warth
##' @export
qarp.pvalue <- function(mat, aVec, permutations=NULL) {
    # do some sanity testing of the arguments.
    if (!is(mat, "matrix")) stop("Must pass a matrix containinf=g read depths across a gene, with transcript position down the rows and each experiments across columns")
    if (!is.factor(aVec)) stop("aVec parameter must be a vector of factors.  It is used to separate the columns of 'mat' into groups.")
    if (length(levels(aVec)) < 2) stop("aVec parameter must have at least two factor levels.")
    if (length(aVec) != ncol(mat))
      warning(sprintf(paste0("Length of aVec (%d) does not match #rows of mat (%d).\n",
                             "Make sure 'mat' has transcript positions in rows and experiments in columns."),
                      length(aVec), nrow(mat)))
    
    distance <- qarp.distance(t(mat))   
 
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
    # available through adonis() then you will want to precompute your own
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
##' At present, this routine is a thin veneer on the 'dist' function.
##' The indirection is useful only because we may want to replace the
##' distance function at some point so it will be conventient if all
##' distance calculations go through one function.
##' @param mat matrix of points, each row is a single point and each
##' column is a single dimensional component of that point.
##' @param method the distance measure to be used. This must be one of
##' "euclidean", "maximum", "manhattan", "canberra", "binary" or
##' "minkowski". Any unambiguous substring can be given.
##' 
##' @return a distance matrix of class 'dist' like this returned by base::dist or vegan::vegdist 
##' @author Chris Warth
##' @export
qarp.distance <- function(mat, ...) {
    distance <- dist(mat, ...)
}


##' Read a tab-separated-value file and cache the resulting dataset so
##' reading it will be faster next time.
##'
##' This routine is meant to work in conjunction with the script
##' \code{bam2depth.py} to access read-depth measurements across
##' genes.  The input to this routine is expected to be a
##' tab-separated text file with the following mandatory columns:
##' 
##' "symbol" - a text field containing the name of a gene or feature.
##' Usually the common gene name is used, but RefSeq or uniprot names
##' could be sued as well.  The same value should be used for any row
##' that contains dep-depth for the same feature.
##'
##' "transPos" - transcript position.  A positive integer indicating
##' position on the transcript where the read depth was measured
##' relative to the start of the feature.  For genes, the firt
##' postition should be 1.
##'
##' "exonNumber" - A positive integer indicating the current exon with
##' the genetic feature.  This column is required but may be set to
##' any single value if exonic information is ot available.
##'
##' Beyond these mandatory columns should be the actual depth data for individual experiments, e.g.
##' symbol	transPos	exonNumber	Benign_1312_501369.depth	Benign_1675_224804.depth
##'
##' These are the columns that must be present (in any order) for the
##' software to function, but it is a non-exclusive list, and other
##' columns besides those dedicated solely to read-depth data may also
##' be present but will be ignored.
##' 
##' @param the name of the file which the data are to be read from.   See \code{help('read.csv')} for more information.
##' @param sep the field separator character.  Default is a single tab character.  See \code{help('read.csv')} for more information.
##' @param cache if TRUE, attempt to read and write a cached version of the file in the same
##' directory; otherwise no cache is used.
##' @return
##' dataframe as parsed from the text file.
##' @author Chris Warth
##' @export
qarp.read.cachedtsv <- function(filename, sep='\t', cache=TRUE) {
    data <- NULL
    base <- sub("^(.*)\\..*", "\\1", filename)
    rdafile = paste0(base, ".Rda")
    options(warn=-1)
    # try to read the ached version
    if (cache) {
        e <- try({
                  supressWarnings(load(paste0(base, ".Rda")))
              }, silent=TRUE)
        if (class(e) != "try-error") {
            attr(data, "filename") = filename
            attr(data, "cachefile") = rdafile
            attr(data, "readcache") = TRUE
        }
    }
    if (is.null(data)) {
        data <- read.csv(filename, sep=sep, stringsAsFactors=FALSE)
        ## HACK - remove the first two rows of totals from the
        ## dataset....  in some datasets these are total number of
        ## reads or total number of alignments.  need to make this more
        ## discriminating and only remove these rows if they in fact
        ## looks like the totals described.
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
##' @return a data frame containing pvalues indexed by gene names
##' @author Chris Warth
##' @export
qarp <- function(data, aVec, genelist=NULL, permutations=NULL) {
    if (is.null(genelist)) {
        ## extract the list of genes.
        genelist = levels(data$symbol)
        stopifnot(!is.null(genelist))
    }
    
    pvalues = sapply(X = genelist,
        FUN = function(gene, data, aVec, permutations) {
            mat <- qarp.matrix(data, gene, aVec)
            qarp.pvalue(mat, aVec, permutations=permutations)
        }, data = data, aVec = aVec, permutations = permutations)
    data.frame(row.names=genelist, pvalues)
}

##' Create a matrix of gene specific read depths from a larger
##' dataset.  The resulting matrix is suitable for use in other QARP
##' functions.
##'
##' The resulting matrix will be restricted to data for a single gene.  
##' 
##' @param data a data frame containing raw read counts for multiple positions along many genes.
##' @param gene the common name of a gene found in the data.
##' @param aVec a list of factors corresponsing to datasets in the
##' data frame.  This is used to identify which columns in the matrix
##' correspond to treated and control datasets.  There should be only
##' two factor levels in the list and each entry should have a name
##' corresponding to the name of one column in the data frame.
##' @return Returns a matrix suitable for use in the other QARP functions like qarp.plotProfile.
##' @author Chris Warth
##' @export
qarp.matrix <- function(data, gene, aVec) {
    stopifnot(all(c("symbol", "transPos", "exonNumber") %in% names(data)))

    if (!gene %in% data$symbol) {
        warning(paste0("feature ", gene," does not exist in symbol column of data frame."))
        return()
    }
    data.gene <- data[data$symbol==gene,]
    data.gene <- data.gene[order(data.gene$transPos),]		# order by transcript position
    rownames(data.gene) = data.gene$transPos
    isplices <- which(diff(data.gene$exonNumber) != 0)
    splices <- data.gene[isplices, 'transPos']
    data.gene = data.gene[, names(aVec), drop=FALSE]
    mat = as.matrix(t(data.gene))

    # if any of the entries end up as NA, simply set them to 0
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


##' plot read profiles across a gene.
##'
##' This routine creates a plot on the current graphics device that
##' depicts a collection of RNA-SEQ profiles drawn over a gene
##' transcript.
##'
##' The profiles will be grouped according to the factors in aVec.
##'
##' If \code{showplices} is true, the spice junctions will be indicated on the plot.
##' These are derived from the read-depth data and may not be accurate for sparsly populated genes.
##' 
##' If \code{invert} is true, one of the groups of profiles will be inverted
##' so it plots below the horizontal axis.   Someone find this a more appealing presentation.
##' 
##' \code{xlim} allows you to specify the limits of
##' the x-axis.
##' @param mat a matrix of read depths for a single gene. Each row is a single nucleotide position on the gene, and each column is seperate experiment.
##' The order of columns should correspond to the order of factors in \code{aVec}.
##' @param aVec list of factors used to group the profiles. 
##' @param showsplices Show splice junctions.
##' @param invert Draw one of the groups of profiles below the horizontal axis.
##' @param xlim optional 2-element vector giving limitations of x-axis.   Using this paramter you can zoom in on particular region of interest.
##' @param ... optional additional parameters to pass to plotting routines (see \code{\link{par}})
##' @author Chris Warth
##' @export
qarp.plotProfile <- function(mat, aVec, showsplices=TRUE, invert=FALSE, ...) {
    # do some sanity testing of the arguments.
    if (!is.factor(aVec)) stop("aVec parameter must be a vector of factors.  It is used to separate the columns of 'mat' into groups.")
    if (length(levels(aVec)) < 2) stop("aVec parameter must have at least two factor levels.")
    if (length(aVec) != ncol(mat))
      warning(sprintf(paste0("Length of aVec (%d) does not match #rows of mat (%d).\n",
                             "Make sure 'mat' has transcript positions in rows and experiments in columns."),
                      length(aVec), nrow(mat)))

    gene <- attr(mat, "gene")
    mypalette <- c("red", "blue")
    m = max(mat)
    if (invert) {
        ylim=c(-m,m)
    } else {
        ylim=c(min(mat),m)
    }

    matplot(mat[,aVec==(levels(aVec)[1])], col=mypalette[1], lty=1, type="l", ylim=ylim, ...)
    if (any(aVec==(levels(aVec)[2])))
      matlines((if (invert) -1 else 1)*mat[,aVec==(levels(aVec)[2])], col=mypalette[2], lty=1, ylim=ylim, ...)
    if (showsplices & !is.null(attr(mat,"junctions"))) {
        j <- attr(mat, "junctions")
        abline(v=j, col=addTrans("red",100), lty=2)
    }

    xpd <- par("xpd")
    par(xpd=TRUE)  # don't clip to the plot area
    smartlegend(x="center", y="top",
            levels(aVec)[1:2], pch=NA, col=mypalette, lty=1, horiz=T, bty="n", inset=1.05)
    par(xpd=xpd)  # reset clipping
    return()
}


##' plot multi-dimensional sampling points.
##'
##' This routine creates a plot on the current graphics device that
##' depicts the distribution of points after multi-dimentional
##' sampling has reduced the dimensionality of the points.  The input
##' matrix contains multiple points, one per row, in a high
##' dimensional space.  e.g. a single RNA-SEQ profile across a
##' 1000-base gene would be one point in a 1000-dimensional space.
##' Multidimensional scaling reduces these points to two dimensions,
##' while maintaining their relative distance from one another.
##' Points close togther in 1000-dimensional space will be close
##' together in 2-dimensional space, and vice versa.
##' @param mat 
##' @param aVec 
##' @author Chris Warth
##' @export
qarp.plotMDS <- function(mat, aVec, ...) {

    distance <- qarp.distance(t(mat))

    # get the two dimension points to be drawn on the scatter plot
    fit <- cmdscale(distance, eig=TRUE, k=2)
    df <- data.frame(fit$points)

    palette <- c("red", "blue")

    # make a factor colum that indicates the condition (treated or
    # control) of each dataset
    cond <- aVec

    # add a column indicating what color should be used for each point.
    df <- transform(df, col=palette[cond])

    #par(mar=c(3, 3, 0.5, 1))  # Trim margin around plot [bottom, left, top, right]

    #par(mgp=c(1.5, 0.2, 0))  # Set margin lines; default c(3, 1, 0) [title,labels,line]
    #par(xaxs="r", yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

    defaults <- list(xaxs = 'r', yaxs = 'r', mgp = c(1.5, 0.2, 0), pch = 21, cex=1.5, yaxt="n")
    args <- modifyList(defaults, list(x =df$X1, y=df$X2, col=as.character(df$col), bg=as.character(df$col), ...))

    do.call("plot", args)

    ## plot(X2 ~ X1, data=df, col=as.character(col), pch=21, bg=as.character(col),
    ##         cex=1.5, frame.plot=F, yaxt="n", ...)

    ## smartlegend(x="center", y="top",
    ##         levels(aVec)[1:2], pch=21, col=mypalette, horiz=T, bty="n", xpd=TRUE, inset=-.1)

    ## ticks = pretty(c(min(df$X2), max(df$X2)), 3)

    ## # par(mgp=c(axis.title.position, axis.label.position, axis.line.position))
    ## mgp <- par("mgp")
    ## mgp[2] <- 0.5
    ## par(mgp=mgp)

    ## axis(2, at=ticks, labels=T, lwd=0, lwd.ticks=1, lty="solid", las=1, cex.axis=0.9)
    ## grid()

    invisible()
}

