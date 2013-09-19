require(vegan)


mds <- function(distance, aVec, permutations=10000) {
    ## How likely is the grouping that we see?   If the datasets were partitioned into
    ## groups in a  random fashion, how often would we see a grouping like the observed grouping?

    ## adonis(formula, data, permutations = 999, method = "bray",
    ##        strata = NULL, contr.unordered = "contr.sum",
    ##        contr.ordered = "contr.poly", ...)

    ## There are two ways to call vegan::adonis.  The first used raw data
    ## and the second uses a distance matrix calculated from the raw data.
    ## In the former case, you pass a matrix of raw measurments to
    ## adonis() and it calculates the distance matrix for you using the
    ## algorithm specified by the 'method' argument.  In the second type
    ## of usage, you precompute a distance matrix from the raw data using
    ## any algorithm you like and then pass the distance matrix as the
    ## first argument to adonis().  adonis() can tell which type of usage
    ## you expect by looking at the type of the first argument; if it is
    ## of class 'dist' then you have precomputed a distance matrix,
    ## otherwise you want adonis() to do it for you.
    ##
    ## Why would you want to use adonis() one way versus the other?  It is
    ## usually going to be easier to pass raw data to adonis() and let it
    ## calculate the distance matrix for you.  However, if you want to try
    ## out a distance calculation that is not one of the algorithms
    ## available through adonis() then you will want to precopute your own
    ## distance matrix and pass it in.

    df <-  adonis(distance ~ aVec, permutations=permutations)
    pvalue <- df$aov.tab[6][1,1]
}

pmds <- function(gene, data = NULL) {
    ## Make a data matrix, where each row is an experiment, and each
    ## column is a transcript position.  The value at each position
    ## represents the read depth at that particular position (col) in one
    ## experiment (row).
    data <- data[data$symbol==gene,]
    data <- data[order(data$transPos),]		# order by transcript position
    rownames(data) = data$transPos
    data = subset(data,,c(-symbol, -transPos))
    mat = as.matrix(t(data))

    ## Calculate a distance matrix and perform multi-dimentsional
    ## sampling to reduce the data to 2-D.  The default distance measure
    ## used by dist() is euclidian, but other diatance measures can also
    ## be used.
    distance <- dist(mat, method="euclid")
    fit <- cmdscale(distance, eig=TRUE, k=2)

    message(paste0("calculating ", gene), appendLF=FALSE)
    aVec <- as.factor(sapply(rownames(mat) %in% control, function(x) if (x) {2} else {3}))
    p <- mds(distance, aVec)
    message(sprintf("\t%.5f",p))
    return(p)
}

filename = "PEO1-RPT-corrUp-coverage.tsv"

## specifiy which datasets are grouped together
control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
datasets = c(control = control, treated = treated)

## read the dataset from a file
message("reading dataset...", appendLF = FALSE)
data = read.csv(filename, sep='\t')
message("done.")
data = data[, c("symbol", "transPos", datasets)]

## extract the list of genes.
genes = levels(data$symbol)

pvalues = sapply(X = genes, FUN = pmds, data = data)




