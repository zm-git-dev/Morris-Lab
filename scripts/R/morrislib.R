suppressMessages( library(RMySQL) )
suppressMessages( library(calibrate) )
suppressMessages( library(affy) )
library(preprocessCore)

suppressMessages( require(maptools) )
suppressMessages( library(plotrix) )


## fetchData - retrieve a dataframe containing the data we are interested in comparing.
## in this case we are retrieving the data from a sql database.
##
## datasets - contains the names of multiple datasets, e.g. c("110112_A_MM1", "110112_B_MM1")
##
## mincount - the minimum number of reads that a gene expression must possess before it will be counted.
##
## group - The mysql database group to use.  "morrisdata" assumes the
## database is accessible on the local network.  "remote" assumes the
## database is accessible on localhost:6607.  Usually this means you
## are at a remote location and you have set up a tunnel via ssh to
## the database at UW.
morris.fetchData <- function(datasets, group=NULL) {

    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
            # con <- dbConnect(drv, user="readonly", password="readonly", dbname="morris", host="localhost")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop(paste0("Could not connect to database: ", e$message));
        }

        getGenome <- function(dataset) {
            ## SQL query for retrieving the genome of each experiment.
            query <- paste0("select genome as '", dataset,"' from datasets_tbl d ",
                     "where d.name like '", dataset, "%'")
            df <- dbGetQuery(con, query)
            if (nrow(df) == 0) {
                print(paste0("Error during database query: no dataset named '", dataset,"'"))
                print("Available datasets:")
                print(morris.datasets(group))
                stop("No Data")
            }
            return(df)
        }

        ## retrieve the genome name for each experiment listed in datasets.
        gs <- lapply(datasets, getGenome)
        gf = Reduce(function(...) merge(..., all=T, suffixes=c("", ".2")), gs)

        ## check to see that each experiment has the same genome!
        ## ( it won't be useful to compare datasets aligned with different genomes. )
        ##
        if (nrow(unique(t(gf[1,]))) != 1) {
            print(paste0("Error during database query: not all selected datasets align to the same genome"))
            print("Selected genomes:")
            print(gs)
            stop("No Data")
        }

        ## If we get here we know all datasets used the same genome so just pick the first one.
        ## This genome will be added as an attribute of the dataframe returns from this function.
        ## (useful when retrieving a set of knowngenes to go with this dataframe.)
        genome=gf[1,1]
        
        f <- function(dataset) {
            ## SQL query for retrieving the per-gene alignment count of a single dataset.
            ## These individual tables will be joined together based on the gene name
            ## as a common key.   Oddly eough it is faster to do this in R than it is
            ## in SQL.
            query <- paste0("select feature as name,count(feature) as '", dataset, "' from ",
                            "new_alignments_tbl a join datasets_tbl d on a.dataset_id=d.id ",
                            "where d.name like '", dataset, "%' group by feature")
            df <- dbGetQuery(con, query)
            if (nrow(df) == 0) {
                print(paste0("Error during database query: no data in dataset '", dataset,"'"))
                print("Available datasets:")
                print(morris.datasets(group))
                stop("No Data")
            }
            return(df)
        }

        ## retrieve the experimental data for each experiment listed in datasets.
        ds <- lapply(datasets, f)

        ## merge all the experimental results into a single dataframe based on
        ## the refseq gene name found in the "refseq" column.
        ##
        ## For details, see
        ## http://stackoverflow.com/questions/8091303/merge-multiple-data-frames-in-a-list-simultaneously
        ## 
        df = Reduce(function(...) merge(...,all=T, by="name", suffixes=c("", ".2")), ds)

        attr(df, "genome") <- genome
        
        # assert that there are no duplicate gene names in the list.
        stopifnot(nrow(df[duplicated(df[,1]),]) == 0)

    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })

    return (df)
}

morris.commonnames <- function(refseq, genome, group=NULL) {
    ## retrieve the knowngene database so we can give common names to the
    ## refseq names used to label each row in our data.
    ## Refseq names have the advantage of being unique, but they are not memorable.
    ##
    kg <- morris.getknowngenes(genome, group=group)

    ## match the refseq names with the common names, using the refseq names as a key
    tmp = kg$name2[match(refseq, kg$name)]

    ## Some refseq names don't have corresponding common names.
    ## For those, reuse the refeseq names as the common name.
    tmp[is.na(tmp)] = refseq[is.na(tmp)]

    return (tmp)

}

## fetch descriptions of datasets from the mysql dtabase.
##
morris.fetchdesc <- function(datasets, group=NULL) {
    descQuery <- paste0("select d.name,e.description from experiments_tbl e ",
                        "join datasets_tbl d on d.expr_id=e.id where ",
                        paste0("d.name like '", datasets, "%'", collapse=" or "))
    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
            ## con <- dbConnect(drv, user="readonly", password="readonly", dbname="morris", host="localhost")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop(paste0("Could not connect to database: ", e$message));
        }

        df <- dbGetQuery(con, descQuery)
        rownames(df) <- df[,1]
        df <- subset(df, select=-c(name) )
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })
    return(df)
}


##
## datasets - contains the names of two datasets, e.g. c("110112_A_MM1", "110112_B_MM1")
##
## mincount - the minimum number of reads that a gene expression must possess before it will be counted.
##
## group - The mysql database group to use.  "morrisdata" assumes the
## database is accessible on the local network.  "remote" assumes the
## database is accessible on localhost:6607.  Usually this means you
## are at a remote location and you have set up a tunnel via ssh to
## the database at UW.

morris.maplot <- function(datasets, mincount=25, group=NULL, normalization="quantile",
                          genes=NULL, report=5) {
    df = morris.genecounts(datasets, mincount=mincount, normalization=normalization, group=group)

    x <- df[(ncol(df)-1):ncol(df)]

    ## identify the entries with the highest and lowest fold changes
    dexpression <- order(x[,2]-x[,1])[1:report]
    dexpression <- append(dexpression, order(x[,1]-x[,2])[1:report])
  
    results <- data.frame(gene=df[dexpression,1],
                          rawA=df[dexpression,2],
                          rawB=df[dexpression,3],
                          normA=x[dexpression,1],
                          normB=x[dexpression,2],
                          lg2A=log2(x[dexpression,1]),
                          lg2B=log2(x[dexpression,2]),
                          diff=(x[,2]-x[,1])[dexpression])
    rownames(results) <- rownames(df)[dexpression]
    results <- results[order(-results$diff),]

    ## nf <- layout(matrix(c(1,2), 1, 2, byrow <- TRUE), widths=c(2,1), heights=c(1,1))
    descriptions = morris.fetchdesc(datasets, group=group)
        
    ma.plot( rowMeans(log2(x)), log2(x[, 1])-log2(x[, 2]), 
            xlab="Mean",ylab="decreased expression - increased expression", cex=0.7) 
    title(paste0(descriptions[,1]," (",rownames(descriptions),")",
                 collapse="\n vs. "))
    textxy(rowMeans(log2(x))[dexpression],(log2(x[, 1])-log2(x[, 2]))[dexpression],
           df[dexpression,'gene'], dcol="red")
  
  
  return (results)
}


##
## Returns a table of the most over-/under-expressed genes.
##
## datasets - the names of datasets in the MySQL database to use, e.g. c("030713_A", "103112_B")
##
## mincount - filter out genes if any dataset has fewer then mincount raw reads before scaling.
##
## normalization - the type of nomalization to apply.
##
## report - this many over/under expressed genes.
##
## genes - extra genes to highlight, in addition to the top over/under expressed, e.g. c("Pbsn", "Sbp")
##
morris.scatter <- function(datasets, mincount=25, group=NULL,
                           normalization="rpkm",
                           genes=NULL, report=5, logscale=FALSE) {

    df = morris.genecounts(datasets, group=group)

    ## normalize the data
    x = morris.normalize(df, normalization, group)

    ## identify entries that do not have proper read depth
    ## and remove them.
    row.sub = apply(df, 1, function(row) (all(!is.na(row)) && all(row >= mincount)))
    x <- x[row.sub,]
    df <- df[row.sub,]
    
    ## add common names for the refseq genes.
    df$common <- morris.commonnames(rownames(df), attr(df, "genome"), group=group)

    ## rearrange the columns to move the last column (common name) to the second column
    ## this is a little tricky b/c number of cols in df is variable.
    df <- df[c(ncol(df), 1:(ncol(df)-1))]

    ## nf <- layout(matrix(c(1,2), 1, 2, byrow <- TRUE), widths=c(2,1), heights=c(1,1))
    descriptions = morris.fetchdesc(datasets, group=group)
    if (logscale==TRUE) {
        plot( x[,1], x[,2], cex=.5, log="xy",
             panel.first = grid(4,4),
             xlab=paste0(descriptions[datasets[1],1], " (log)"),
             ylab=paste0(descriptions[datasets[2],1], "(log)"))
        abline(lm(log10(x[,2])~0+log10(x[,1])), col="red")

    } else {
        plot( x[,1], x[,2], cex=.5,
             xlab=descriptions[datasets[1],1],
             ylab=descriptions[datasets[2],1])
        abline(lm(x[,2]~0+x[,1]), col="blue")
    }
    title(paste0(descriptions[,1]," (",rownames(descriptions),")",
                 collapse="\n vs. "), cex=.7)

    ## Prepare a dataframe with information about the most differentially
    ## expressed points.
    dexpression <- order(x[,2]-x[,1])[1:report]
    dexpression <- append(dexpression, order(x[,1]-x[,2])[1:report])
    dexpression <- append(dexpression, match(genes, df[,1]))
  
    results <- data.frame(gene=df[dexpression,1],
                          rawA=df[dexpression,2],
                          rawB=df[dexpression,3],
                          normA=x[dexpression,1],
                          normB=x[dexpression,2],
                          diff=(x[,2]-x[,1])[dexpression])
    rownames(results) <- rownames(df[dexpression,])
    results <- results[order(-results$diff),]

    ## label the most differentialy expressed points.
    text(x[dexpression,],  df[dexpression,'common'], adj=c(-0, -0.3), cex=.7, col="red")

    return (results)
}


## A function to use identify to select points, and overplot the
## points with another symbol as they are selected
identifyPch <- function(x, y = NULL, n = length(x), pch = 19, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch)
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    return(res)
}



## retrieve a list of dataset names that have per-gene expression values stored in the database.

## This routine will let you fetch the names of datasets based on
## criteria such as organism, tissue, or genotype.
##
morris.datasets <- function(group=NULL, organism=NULL, tissue=NULL, genotype=NULL) {

    ## a query to find all datasets that have alignments in the new_alignments_tbl table.
    query.alignments <- paste("select d.name, d.expr_id from",
                              "(SELECT dataset_id FROM morris.new_alignments_tbl",
                              "group by dataset_id) a",
                              "join datasets_tbl d ON d.id = a.dataset_id")

    ## a query to find experiments that use a certain organism, tissue, or genotype.
    query.experiments <- "select id from experiments_tbl"
    selects <- c( ifelse(!is.null(organism),  paste0("organism like '", organism, "'"), NA),
                  ifelse(!is.null(tissue),  paste0("tissue like '", tissue, "'"), NA),
                  ifelse(!is.null(genotype),  paste0("genotype like '", genotype, "'"), NA) )
    if (any(!is.na(selects))) 
        query.experiments = paste( query.experiments, "where", paste(selects[!is.na(selects)], collapse=" and "))

    ## Assemble the complete query from the sub-queries.
    query <- paste("select d.name from (", query.alignments, ") d",
                   "join (", query.experiments, ") e",
                   "ON d.expr_id = e.id")
    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        df <- dbGetQuery(con, query)
    }, warning = function(w) {
        print(w$message)
    }, error = function(e) {
        print(e$message)
        traceback(2)
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })

    return (df$name);
}


morris.genecounts <- function(datasets, group=NULL) {
    ds <- morris.fetchData(datasets, group)
    rownames(ds) <- ds[,1]
    df <- subset(ds, select=-c(name) )
    
    ## Copy the attribute for the original data.frame if necessary
    ## see http://stackoverflow.com/a/10420036/1135316
    n = names(df)
    mostattributes(df) = attributes(ds)
    names(df) = n
    
    return(df)
}


## Retrieve the list of known refseq genes from the database.
##
morris.getknowngenes <- function(genome, gene=NULL, group=NULL) {

    # SQL query for retrieving the knowngene list.  The proper genome is selected from
    # the first dataset.
    kquery <- paste0("select * from knowngenes2 k ",
                     "where k.genome=\"", genome, "\"");
    if (!is.null(gene)) {
        kquery <- paste0(kquery, " and ( name in (", paste0("'",gene,"'", collapse=","), ")")
        kquery <- paste0(kquery, " or name2 in (", paste0("'",gene,"'", collapse=","), ")  )")
    }

    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
            # con <- dbConnect(drv, user="readonly", password="readonly", dbname="morris", host="localhost")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop(paste0("Could not connect to database: ", e$message));
        }


        ## retrieve the knowngene database.
        ## Calculate the length of exonic sequence
        ##
        kg <- dbGetQuery(con, kquery)
        kg <- kg[!duplicated(kg$name),]
        fun <- function(x,y) {sum(as.numeric(y)-as.numeric(x))}
        kg$exonLen <- mapply(fun, strsplit(kg$exonStarts,","), strsplit(kg$exonEnds,","))

    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })
    return(kg)
}



morris.normalize <- function(df, 
                             normalization=c("rpm", "scaled", "quantile", "rpkm"),
                             group=NULL) {
    
    ## normalize to total reads in each experiment.
    normalization <- match.arg(normalization)
    df.cols <- sapply(df, is.numeric)  #select numeric variables
    if (normalization == "rpkm") {

        ## FPKM = 10e9 * C / (N * L)
        ## C = # mapped reads that fell into exons of a particular gene
        ## N = total mapped reads for this experiment
        ## L = exonic base pairs in this gene
        ## http://www.biostars.org/p/11378/
        ## http://www.biostars.org/p/6694/
        ##
        kg <- morris.getknowngenes(attr(df, "genome"), group=group)
        N <- colSums(df[df.cols], na.rm = TRUE)
        L <- kg[match(rownames(df),kg$name),"exonLen"]
        C <- df[df.cols]

        ## sapply gives a matrix.
        ## lapply would return a list of lists (not useful here)
        #tmp <- sapply(N,function(n) {return(n*L)})
        tmp = (N*L)
        tmp <- (10^9)/tmp
        x <- C*tmp
        
    } else if (normalization == "rpm") {
        N <- colSums(df[df.cols], na.rm = TRUE)
        C <- df[df.cols]
        
        x <- C/(N/(10^6))
      
    } else if (normalization == "scaled") {
        N <- colSums(df[df.cols], na.rm = TRUE)
        C <- df[df.cols]

        x <- C/N

    } else if (normalization == "quantile") {
        
        x <- data.frame(normalize.quantiles(data.matrix(df[df.cols])))
        names(x) <- datasets
        rownames(x) <- rownames(df)
    
    }
    return(x)
}
