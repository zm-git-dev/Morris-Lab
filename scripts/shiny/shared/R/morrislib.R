suppressMessages( library(RMySQL) )
suppressMessages( library(calibrate) )
suppressMessages( library(affy) )
library(preprocessCore)

##suppressMessages( require(maptools) )
suppressMessages( library(plotrix) )

options(deparse.max.lines=30)
options(max.print=30)
## options(error = quote({dump.frames(to.file=TRUE); q()}))



##' Retrieve per-gene alignment
##' counts for a collection of datasets.
##'
##' This data is useful for differential gene expression analysis, but
##' it not useful for ribosome profile analysis - use getalignments
##' for that.
##' @title morris.fetchdata
##' @param datasets contains the names of multiple datasets,
##' e.g. c("110112_A_MM1", "110112_B_MM1")
##' @param group The mysql database group to use.  "morrisdata"
##' assumes the database is accessible on the local network.  "remote" assumes the
##' database is accessible on localhost:6607.  Usually this means you
##' are at a remote location and you have set up a tunnel via ssh to
##' the database at UW.
##' @return Returns data frame containing one row per gene, and one
##' column per dataset.  Each entry in the row is a raw read count for
##' one gene in one dataset.  Missing values will be marked NA.  A
##' gene must have a non-zero value in at least one dataset in order
##' to be listed.
##' @author chris
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
            stop("Could not connect to database: ", e$message);
        }

        ## retrieve the genome name for each experiment listed in datasets.
        gs <- morris.getGenome(datasets)

        ## check to see that each experiment has the same genome!
        ## ( it won't be useful to compare datasets aligned with different genomes. )
        ##
        if (length(unique(gs)) != 1) {
            stop("not all selected datasets align to the same genome\n",
                 "Selected genomes:\n", as.character(gs))
        }

        ## If we get here we know all datasets used the same genome so just pick the first one.
        ## This genome will be added as an attribute of the dataframe returns from this function.
        ## (useful when retrieving a set of knowngenes to go with this dataframe.)
        genome=unique(gs)
        
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
                stop("no data in dataset '", dataset,"'\n",
                     "Available datasets:\n",
                     as.character(morris.datasets(group)), call.=TRUE)
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

    cn = data.frame(refseq, row.names=1)

    ## match the refseq names with the common names, using the refseq names as a key
    cn$common = kg$name2[match(refseq, kg$name)]

    ## Some refseq names don't have corresponding common names.
    ## For those, reuse the refeseq names as the common name.
    rows.na = is.na(cn$common)
    cn$common[rows.na] = refseq[rows.na]

    return (cn)

}


## Fetch descriptive information about a group of
## datasets.  This routine builds a dataframe from a combination of
## the datasets_tbl and the experiments_tbl.
##
morris.fetchinfo <- function(datasets, group=NULL) {
    query = paste0("SELECT e.name as 'experiment',d.name as 'dataset',",
                   "       e.subjectID,e.organism,e.age,e.description,",
                   "       e.tissue,e.genotype ",
                   "FROM morris.datasets_tbl d ",
                   "JOIN experiments_tbl e ON d.expr_id = e.id where ",
                   paste0("d.name like '", datasets, "%'", collapse=" or "))
    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop("Could not connect to database: ", e$message);
        }

        df <- dbGetQuery(con, query)
        rownames(df) <- df[,'dataset']
        df <- subset(df, select=-c(dataset) )
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })
    return(df)
}



## Fetch statistical information about a group of
## datasets.  This routine builds a dataframe from a combination of
## the datasets_tbl and the experiments_tbl.
##
morris.fetchstats <- function(datasets, group=NULL) {
    query = paste0("SELECT e.name as 'experiment',d.name as 'dataset',",
                   "       e.raw_count, d.trimmed_count,",
                   "       d.aligned_count, d.in_uniq_exons_count,",
                   "       (d.trimmed_count - d.nonrrna_count) as 'rRNA'",
                   "FROM morris.datasets_tbl d ",
                   "JOIN experiments_tbl e ON d.expr_id = e.id where ",
                   paste0("d.name like '", datasets, "%'", collapse=" or "))
    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop("Could not connect to database: ", e$message);
        }

        df <- dbGetQuery(con, query)
        rownames(df) <- df[,'dataset']
        df <- subset(df, select=-c(dataset) )
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })
    return(df)
}



## fetch descriptions of datasets from the mysql dtabase.
## FIXME - this should be replaced with morris.fetchinfo
##
morris.fetchdesc <- function(datasets, group=NULL) {
    .Deprecated("morris.fetchinfo", "morrislab", "Fetchdesc is superceded by fetchinfo.")
    df = morris.fetchinfo(datasets, group=group)
    return(df[, c("description"), drop=FALSE])
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
    descriptions = morris.fetchinfo(datasets)[,"description", drop=FALSE]
        
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
    descriptions = morris.fetchinfo(datasets)[,"description", drop=FALSE]
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
    selects <- c( ifelse(!is.null(organism),  paste0("organism = '", organism, "'"), NA),
                  ifelse(!is.null(tissue),  paste0("tissue like '", tissue, "'"), NA),
                  ifelse(!is.null(genotype),  paste0("genotype = '", genotype, "'"), NA) )
    if (any(!is.na(selects))) 
        query.experiments = paste( query.experiments, "where", paste(selects[!is.na(selects)], collapse=" and "))

    ## Assemble the complete query from the sub-queries.
    query <- paste("select d.name from (", query.alignments, ") d",
                   "join (", query.experiments, ") e",
                   "ON d.expr_id = e.id")

    if (is.null(group)) 
        group="remote"

    result <- tryCatch({
        con <- dbConnect(dbDriver("MySQL"), group=group)
        df <- dbGetQuery(con, query)
    }, warning = function(w) {
        print(w$message)
    }, error = function(e) {
        print(e$message)
        traceback(2)
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
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
    stopifnot(!is.null(genome))
    
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
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop("Could not connect to database: ", e$message);
        }


        ## retrieve the knowngene database.
        ## Calculate the length of exonic sequence
        ##
        kg <- dbGetQuery(con, kquery)   # may return NULL!
        ##stopifnot(nrow(kg) > 0)
        message("back from dbGetQuery")
        if (!is.null(kg) && nrow(kg)) {
            # FIXME - do we reallt want to filter dusplicated names at this low level?
            # probably not!

            ##print(paste("duplicated = ", duplicated(kg$name)))
            kg <- kg[!duplicated(kg$name),]

            # sum up all the exons of these genes to get a total transcript length.
            fun <- function(x,y) {sum(as.numeric(y)-as.numeric(x))}
            kg$exonLen <- mapply(fun, strsplit(kg$exonStarts,","),
                                 strsplit(kg$exonEnds,","))
        }
        message("returning from getknowngenes")
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })
    return(kg)
}



morris.normalize <- function(df, 
                             normalization=c("TC", "rpm", "scaled", "quantile", "rpkm"),
                             genoe=NULL) {
    ## normalize to total reads in each experiment.
    normalization <- match.arg(normalization)
    df.cols <- sapply(df, is.numeric)  #select numeric variables

    if (normalization == "rpkm") {
        if (is.null(genome)) {
            stopifnot(!is.null(attr(df, "genome")))
            genome = attr(df, "genome")
        }

        ## FPKM = 10e9 * C / (N * L)
        ## C = # mapped reads that fell into exons of a particular gene
        ## N = total mapped reads for this experiment
        ## L = exonic base pairs in this gene
        ## http://www.biostars.org/p/11378/
        ## http://www.biostars.org/p/6694/
        ##
        kg <- morris.getknowngenes(genome, group=group)
        N <- colSums(df[df.cols], na.rm = TRUE)
        L <- kg[match(rownames(df),kg$name),"exonLen"]
        C <- df[df.cols]

        ## sapply gives a matrix.
        ## lapply would return a list of lists (not useful here)
        ## as tempting as it might be, do NOT replace this with,
        ## tmp = (N*L)
        ## that does not work if the dataframe has multiple datasets in it.
        tmp <- sapply(N,function(n) {return(n*L)})
        tmp <- (10^9)/tmp
        x <- C*tmp
    } else if (normalization == "rpm") {

        ## Normalize to mapped reads per million mapped reads
        ##  http://stackoverflow.com/questions/13830979/#13831155
        C <- df[df.cols]
        N <- colSums(C, na.rm = TRUE)

        ## x = scale(C, center=FALSE, scale=colSums(C)/(10^6)
        x = sweep(C,2,colSums(C)/(10^6),`/`)
        
    } else if (normalization == "scaled") {
        C <- df[df.cols]
        N <- colSums(C, na.rm = TRUE)

        ## here is an alternative technique....
        ## df = sweep(df,2,(colSums(df)/(10^6)), '/')
        x = sweep(C,2,colSums(C),`/`)

    } else if (normalization == "quantile") {
        
        x <- data.frame(normalize.quantiles(data.matrix(df[df.cols])))
        names(x) <- datasets
        rownames(x) <- rownames(df)
    
    } else if (normalization == "TC") {  
        ## Total Count
        ## Gene counts are divided by the total number of mapped reads
        ## (or library size) associated with their lane and multiplied
        ## by the mean total count across all the samples of the
        ## dataset.
        C <- df[df.cols]
        N <- colSums(C, na.rm = TRUE)

        C = sweep(C,2,N,`/`)
        x = sweep(C,2,mean(N),`*`)
    }
    return(x)
}


##' Retrieve alignments for a single dataset from the database.  If gene is specified, only
##' retrieve alignments for that particular gene.
##'
##' . details
##' @param dataset  A string representation of a dataset name from the SQL table 'datasets_tbl'
##' @param gene A string name of a gene, either common name or RefSeq name
##' @param group Optional specification of the group used to connect to the MySQL database.
##' @return dataframe containing alignments, on alignment per row.
##' @author Chris Warth
morris.getalignments <- function(dataset, gene=NULL, group=NULL) {
    ## this routine operates only on ONE dataset name - NOT a vector!
    stopifnot(length(dataset) == 1)

    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop("Could not connect to database: ", e$message);
        }

        ## check to see that each experiment has the same genome!
        ## it won't be useful to compare datasets aligned with different genomes. 
        ## ( actually this routine only works with a single dataset, so really there is no
        ##   need to consolidate the results from multiple datasets. )
        gs <- morris.getGenome(dataset)
        genome <- unique(gs)
        if (length(genome) == 0) {
            stop("Error during database query: no dataset found with the name ", paste0(dataset, collapse=", "))
        } else if (length(genome) != 1) {
            stop("Error during database query: not all selected datasets align to the same genome")
        }

        ## SQL query for retrieving the per-gene alignment count of a single dataset.
        ## These individual tables will be joined together based on the gene name
        ## as a common key.   Oddly eough it is faster to do this in R than it is
        ## in SQL.
        query <- paste0("select feature,position,length from ",
                        "new_alignments_tbl a join datasets_tbl d on a.dataset_id=d.id ",
                        "where d.name like '", dataset, "%'")
        if (!is.null(gene)) {
            query <- paste0(query, " and a.feature in (", paste0("'",gene,"'", collapse=","), ")")
        }
        df <- dbGetQuery(con, query)
        if (nrow(df) == 0) {
            stop("Error during database query: no alignments in dataset '", dataset,
                 "' for gene(s) ", paste0(gene, collapse=","), "\n",
                 "Available datasets:\n",
                 as.character(morris.datasets(group)), call.=TRUE)
        }

        attr(df, "genome") <- genome
        attr(df, "dataset") <- dataset
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })

    return(df)
}


## Retrieve the genome used to align a setset (or datasets)
## This routine take a list of datasets or a single dataset string
## and returns a list of strings (or single string) representing the
## genome used to align this dataset.
## e.g. morris.getGenome("110313_A") => "hg18"
##      morris.getGenome(c("110313_A", "103113_B")) => ("hg18", "mm9")
##
morris.getGenome <- function(dataset, group=NULL) {
    df <- NULL
    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop("Could not connect to database: ", e$message);
        }

        ## SQL query for retrieving the genome of each experiment.
        query <- paste0("select name,genome from datasets_tbl d ",
                        "where ",
                        paste0("name like '", dataset, "%'", collapse=" or "))
        df <- dbGetQuery(con, query)
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })
    x = pmatch(dataset, df[,1])
    if (ncol(df))
      return(df[x,2])
    else
      return(NULL)
}
