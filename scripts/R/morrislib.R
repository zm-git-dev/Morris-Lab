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
    # SQL query for retrieving the knowngene list.  The proper genome is selected from
    # the first dataset.
    kquery <- paste0("select distinct name,name2 as gene from knowngenes2 k ",
                     "join (select genome from experiments_tbl e join datasets_tbl d ",
                     "on d.expr_id=e.id where d.name like '", datasets[1],
                     "%') g on g.genome=k.genome");

    # SQL query for retrieving the set of datasets with alignments.
    dsquery <- paste0("select d.name from (SELECT dataset_id FROM morris.new_alignments_tbl ",
                      "group by dataset_id) a join datasets_tbl d on d.id=a.dataset_id;")

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

        f <- function(dataset) {
            # SQL query for retrieving the the alignment data of a single dataset.
            query <- paste0("select feature as name,count(feature) as '", dataset, "' from ",
                            "new_alignments_tbl a join datasets_tbl d on a.dataset_id=d.id ",
                            "where d.name like '", dataset, "%' group by feature")
            df <- dbGetQuery(con, query)
            if (nrow(df) == 0) {
                print(paste0("Error during database query: no data in dataset '", dataset,"'"))
                print("Available datasets:")
                print(morris.datasets(group));
                stop("No Data");
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

        # assert that there are no duplicate gene names in the list.
        stopifnot(nrow(df[duplicated(df[,1]),]) == 0)

        ## retrieve the knowngene database so we can give common names to the
        ## refseq names used to label each row in our data.
        ## Refseq names have the advantage of being unique, but they are not memorable.
        ##
        knowngenes <- dbGetQuery(con, kquery)

        # join the knowngene list with our sample data using the refseq gene name column.
        # keep just the refseq name, the cmmon gene name, and the two columns of count data.
        df <- merge(df,knowngenes, all.x=TRUE, by="name")
        
        ## rearrange the columns to move the last column (common name) to the second column
        ## this is a little tricky b/c number of cols in df is variable.
        df = df[c(1, ncol(df),2:(ncol(df)-1))]

        # Some refseq names don't have corresponding common names.
        # For those, reuse the refeseq names as the common name.
        df[is.na(df[,2]),2] <- df[is.na(df[,2]),1]
        
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })

    return (df)
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
    df = morris.genecounts(datasets, mincount=mincount, group=group)

    ## normalize to total reads in each experiment.
    normalization <- match.arg(normalization)
    x<- switch(normalization,
               scaled = df[2:ncol(df)]/colSums(df[2:ncol(df)]),
               quantile = data.frame(normalize.quantiles(data.matrix(df[2:ncol(df)]))))

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
## normlization - the type of nomalization to apply.
##
## report - this many over/under expressed genes.
##
## genes - extra genes to highlight, in addition to the top over/under expressed, e.g. c("Pbsn", "Sbp")
##
morris.scatter <- function(datasets, mincount=25, group=NULL,
                           normalization=c("rpm", "scaled", "quantile"),
                           genes=NULL, report=5, logscale=FALSE) {
    df = morris.genecounts(datasets, mincount=mincount, group=group)

    ## normalize to total reads in each experiment.
    normalization <- match.arg(normalization)
    x<- switch(normalization,
               rpm = df[2:ncol(df)]/(colSums(df[2:ncol(df)]/(10^6))),
               scaled = df[2:ncol(df)]/colSums(df[2:ncol(df)]),
               quantile = data.frame(normalize.quantiles(data.matrix(df[2:ncol(df)]))))

    ## identify the entries with the highest and lowest fold changes
    dexpression <- order(x[,2]-x[,1])[1:report]
    dexpression <- append(dexpression, order(x[,1]-x[,2])[1:report])
    dexpression <- append(dexpression, match(genes, df[,1]))
    
  
    results <- data.frame(gene=df[dexpression,1],
                          rawA=df[dexpression,2],
                          rawB=df[dexpression,3],
                          normA=x[dexpression,1],
                          normB=x[dexpression,2],
                          diff=(x[,2]-x[,1])[dexpression])
    rownames(results) <- rownames(df)[dexpression]
    results <- results[order(-results$diff),]

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
        abline(lm(x[,2]~0+x[,1]), col="red")
    }
    title(paste0(descriptions[,1]," (",rownames(descriptions),")",
                 collapse="\n vs. "))
    text(x[dexpression,],  df[dexpression,'gene'], adj=c(-0, -0.3), cex=.7, col="red")

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
##
morris.datasets <- function(group=NULL) {
    query <- paste0("select d.name from (SELECT dataset_id FROM morris.new_alignments_tbl ",
                   "group by dataset_id) a join datasets_tbl d on d.id=a.dataset_id;")
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

    return (df);
}


morris.genecounts <- function(datasets, mincount=NA, group=NULL,
                           normalization= c( "scaled", "quantile")) {
    df <- morris.fetchData(datasets, group)
    rownames(df) <- df[,1]
    df <- subset(df, select=-c(name) )

    if (!is.na(mincount)) {
        ## Filter NA values seperately because comparison operators don't work on NA.
        ## (NA > anything is always NA).
        df <- df[complete.cases(df[,2:ncol(df)]),]
        
        ## Filter out rows that contain any raw read count less than 'mincount' for any experiment.
        df = df[rowSums(df[,2:ncol(df)] >= mincount) == ncol(df)-1,]
    }
    return(df)
}
