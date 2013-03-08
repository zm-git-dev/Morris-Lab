suppressMessages( library(RMySQL) )
suppressMessages( library(calibrate) )
suppressMessages( library(affy) )
library(preprocessCore)

suppressMessages( require(maptools) )
suppressMessages( library(plotrix) )


## fetchData - retrieve a dataframe containing the data we are interested in comparing.
## in this case we are retrieving the data from a sql database.
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
morris.fetchData <- function(datasets, group=NULL) {
    # SQL query for retrieving the knowngene list.  The proper genome is selected from
    # the first dataset.
    kquery <- paste0("select distinct name,name2 from knowngenes2 k ",
                     "join (select genome from experiments_tbl e join datasets_tbl d ",
                     "on d.expr_id=e.id where d.name like '", datasets[1],
                     "%') g on g.genome=k.genome");

    # SQL query for retrieving the set of datasets with alignments.
    dsquery <- paste0("select d.name from (SELECT dataset_id FROM morris.new_alignments_tbl ",
                      "group by dataset_id) a join datasets_tbl d on d.id=a.dataset_id;")

    descQuery <- paste0("select d.name,e.description from experiments_tbl e ",
                        "join datasets_tbl d on d.expr_id=e.id where ",
                        paste0("d.name like '", datasets, "%'", collapse=" or "))
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
            query <- paste0("select feature as name,count(feature) as count from ",
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

        ds <- lapply(datasets, f)
        df <- merge(ds[[1]], ds[[2]], by="name", all=TRUE)

        # assert that there are no duplicate gene names in the list.
        stopifnot(nrows(df[duplicated(df[,1]),]) == 0)

        knowngenes <- dbGetQuery(con, kquery) 
        df.data <- merge(df,knowngenes, all.x=TRUE, by="name")[c(1,4,2,3)]
    
        df.desc <- dbGetQuery(con, descQuery)
        rownames(df.desc) <- df.desc[,1]
        df.desc <- subset(df.desc, select=-c(name) )
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })

    return ( list(df=df.data, descriptions=df.desc) )
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

morris.maplot <- function(datasets, mincount=25, group=NULL, normalization="quantile") {
  datalist <- morris.fetchData(datasets, group)
  df <- datalist$df
  
  ## identify each row with the gene name.
  rownames(df) <- df[,1]
  genenames <- df[,2]
  df <- subset(df, select=-c(name, name2) )

  # filter out rows that contain NA.
  df <- df[complete.cases(df), ]
  
  # Filter out rows that contain any raw read count less than 'mincount' for ay experiment.
  df <- df[rowSums(df[,1:2] >= mincount)==2, ]
  
  
  ## normalize to total reads in each experiment.
  if (normalization=="quantile") {
      x <- normalize.quantiles(data.matrix(df[complete.cases(df),]))
  } else if (normalization == "scale") {
      x <- df/colSums(df) 
  } else {
      stop("Unrecognized 'normalization' value: ", normalization)
  }
  
  ## identify the entries with the highest and lowest fold changes
  dexpression <- order(x[,2]-x[,1])[1:5]
  dexpression <- append(dexpression, order(x[,1]-x[,2])[1:5])
  
  results <- data.frame(gene=genenames[dexpression],
                       rawA=df[dexpression,1],
                       rawB=df[dexpression,2],
                       normA=x[dexpression,1],
                       normB=x[dexpression,2],
                       lg2A=log2(x[dexpression,1]),
                       lg2B=log2(x[dexpression,2]),
                       diff=(x[,2]-x[,1])[dexpression])
  
  
  ## make a copy of the results table ti display on the graph.
  display <- results
  display$normA <- format(display$normA, digits=2, nsmall=2)
  display$normB <- format(display$normB, digits=2, nsmall=2)
  display$lg2A <- format(display$lg2A, digits=2, nsmall=2)
  display$lg2B <- format(display$lg2B, digits=2, nsmall=2)
  display$diff <- format(display$diff, digits=2, nsmall=2)
  
  ## nf <- layout(matrix(c(1,2), 1, 2, byrow <- TRUE), widths=c(2,1), heights=c(1,1))

  ma.plot( rowMeans(log2(x)), log2(x[, 1])-log2(x[, 2]), 
           xlab="Mean",ylab=paste0(datasets, collapse=" - "), cex=0.7) 
  title(paste0(datalist$descriptions$description," (",rownames(datalist$descriptions),")",
               collapse="\n vs. "))
  textxy(rowMeans(log2(x))[dexpression],(log2(x[, 1])-log2(x[, 2]))[dexpression],genenames[dexpression])
  
  
  return (results)
}



morris.scatter <- function(datasets, mincount=25, group=NULL, normalization="quantile") {
  datalist <- morris.fetchData(datasets, group)
  df <- datalist$df
  
  ## identify each row with the gene name.
  rownames(df) <- df[,1]
  genenames <- df[,2]
  df <- subset(df, select=-c(name, name2) )

  # filter out rows that contain NA.
  df <- df[complete.cases(df),]
  
  # Filter out rows that contain any raw read count less than 'mincount' for ay experiment.
  df <- df[rowSums(df[,1:2] >= mincount)==2,]
  
  
  ## normalize to total reads in each experiment.
  if (normalization=="quantile") {
      x <- normalize.quantiles(data.matrix(df[complete.cases(df),]))
  } else if (normalization == "scale") {
      x <- df/colSums(df) 
  } else {
      stop("Unrecognized 'normalization' value: ", normalization)
  }
  
  ## identify the entries with the highest and lowest fold changes
  dexpression <- order(x[,2]-x[,1])[1:5]
  dexpression <- append(dexpression, order(x[,1]-x[,2])[1:5])
  
  results <- data.frame(gene=genenames[dexpression],
                       rawA=df[dexpression,1],
                       rawB=df[dexpression,2],
                       normA=x[dexpression,1],
                       normB=x[dexpression,2],
                       lg2A=log2(x[dexpression,1]),
                       lg2B=log2(x[dexpression,2]),
                       diff=(x[,2]-x[,1])[dexpression])
  
  
  ## make a copy of the results table ti display on the graph.
  display <- results
  display$normA <- format(display$normA, digits=2, nsmall=2)
  display$normB <- format(display$normB, digits=2, nsmall=2)
  display$lg2A <- format(display$lg2A, digits=2, nsmall=2)
  display$lg2B <- format(display$lg2B, digits=2, nsmall=2)
  display$diff <- format(display$diff, digits=2, nsmall=2)
  
  ## nf <- layout(matrix(c(1,2), 1, 2, byrow <- TRUE), widths=c(2,1), heights=c(1,1))

  plot( x[,1], x[,2],
       xlab=datalist$descriptions$description[0],
       ylab=datalist$descriptions$description[1])
  abline(lm(x[,2]~x[,1]), col="red");
  title(paste0(datalist$descriptions$description," (",rownames(datalist$descriptions),")",
               collapse="\n vs. "))
  textxy(x[dexpression,1],x[dexpression,2],genenames[dexpression])
  
  
  return (results)
}

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
