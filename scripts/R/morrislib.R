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
morris.fetchData <- function(datasets, mincount, group=NULL) {
    query=paste0("select feature as name,count(feature) as count from ",
      "new_alignments_tbl a join datasets_tbl d on a.dataset_id=d.id ",
      "where d.name like '", datasets, "%' group by feature")

    drv = dbDriver("MySQL")
    if (is.null(group)) {
        message("group is NULL, using default parameters to connect to database...")
        con = dbConnect(drv, group="remote")
#         con = dbConnect(drv, user="readonly", password="readonly",
#           dbname="morris", host="localhost")
    } else  {
        message("group is  NOT NULL...")
        con = dbConnect(drv, group=group)
    }
    f = function(q) { dbGetQuery(con, q) }
    ds = lapply(query, f)
    df = merge(ds[[1]], ds[[2]], by="name", all=TRUE)
    head(df[duplicated(df[,1]),])
    kquery = paste0("select distinct name,name2 from knowngenes2 k ",
      "join (select genome from experiments_tbl e join datasets_tbl d ",
      "on d.expr_id=e.id where d.name like '", datasets[1],
      "%') g on g.genome=k.genome");
    knowngenes = dbGetQuery(con, kquery) 
    df.data = merge(df,knowngenes, all.x=TRUE, by="name")[c(1,4,2,3)]

    
    

    descQuery <- paste0("select d.name,e.description from experiments_tbl e ",
                        "join datasets_tbl d on d.expr_id=e.id where ",
                        paste0("d.name like '", datasets, "%'", collapse=" or "))
    df.desc <- dbGetQuery(con, descQuery)
    print(df.desc)
    rownames(df.desc) = df.desc[,1]
    print(df.desc)
    df.desc <- subset(df.desc, select = -c(name) )
    
    return ( list(df=df.data, descriptions=df.desc) )
}
    

morris.oldfetchData <- function(datasets, mincount, group=NULL) {

    ## Build an SQL query string for retrieving alignment data from the database.
    ##
    # Do an erzats FULL JOIN. of tbl1 and tbl2
    # Because MySQL does not implement full joins, we have to UNION a LEFT JOIN with a RIGHT JOIN.


    # Need to join this result with the knowngenes2 table....csw 2/19/2013

    dataQuery <- paste0('select  COALESCE(k.name2,a.name) as genename, a.* from (SELECT t1.feature as name, ',
                        paste0('t', seq(datasets), '.count as c', seq(datasets),collapse=", "), " from\n",
                        paste0("(select feature,count(feature) as count from new_alignments_tbl a\n",
                               " join datasets_tbl d on a.dataset_id=d.id\n",
                               " where d.name like '", datasets,"%' group by feature) t",
                               seq(datasets), collapse="\nleft outer join \n"),
                        "\non ", paste0("t", seq(datasets), ".feature", collapse="="),
                        "\n UNION \n",
                        "SELECT t2.feature as name, ",
                        paste0('t', seq(datasets), '.count as c', seq(datasets),collapse=","), " from\n",
                        paste0("(select feature,count(feature) as count from new_alignments_tbl a\n",
                               " join datasets_tbl d on a.dataset_id=d.id\n",
                               " where d.name like '", datasets,"%' group by feature) t",
                               seq(datasets), collapse="\nright outer join \n"),
                        " on ", paste0("t", seq(datasets), ".feature", collapse="="),
                        ") a  left outer join knowngenes2 k on k.name=a.name \n")
    
    ## Build an SQL query string for retrieving the descriptions of
    ## experiments from the database.  These are usually used to add
    ## useful labels to axes on plots.
    descQuery <- paste0("select d.name,e.description from experiments_tbl e join datasets_tbl d on d.expr_id=e.id where ",
                    paste0("d.name like '", datasets, "%'", collapse=" or "))

    result = tryCatch({
        drv = dbDriver("MySQL")
        if (is.null(group)) {
            message("group is NULL, using default parameters to connect to database...")
            con = dbConnect(drv, user="readonly", password="readonly",
              		    dbname="morris", host="localhost")
        } else  {
            message("group is  NOT NULL...")
            con = dbConnect(drv, group=group)
        }
        cat(dataQuery)
        df.data = dbGetQuery(con, dataQuery)
        print(head(df.data))
        
        df.desc <- dbGetQuery(con, descQuery)
        print(df.desc)
        rownames(df.desc) = df.desc[,1]
        print(df.desc)
        df.desc <- subset(df.desc, select = -c(name) )
    }, warning = function(w) {
        print(w$message)
    }, error = function(e) {
        print(e$message)
        traceback(2)
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
       ## if (exists("drv")) 
         ## dbUnloadDriver(drv)
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
  datalist <- morris.fetchData(datasets, mincount, group)
  df = datalist$df
  df.desc = datalist$descriptions
  
  ## identify each row with the gene name.
  rownames(df) = df[,1]
  genenames = df[,2]
  df <- subset(df, select = -c(name, name2) )
  
  ## normalize to total reads in each experiment.
  if (normalization=="quantile") {
      x <- normalize.quantiles(data.matrix(df[complete.cases(df),]))
  } else if (normalization == "scale") {
      x <- df/colSums(df) 
  } else {
      stop("Unrecognized 'normalization' value: ", normalization)
  }
  
  ## identify the entries with the highest and lowest fold changes
  dexpression = order(x[,2]-x[,1])[1:5]
  dexpression = append(dexpression, order(x[,1]-x[,2])[1:5])
  
  results = data.frame(gene=genenames[dexpression],
                       rawA=df[dexpression,1],
                       rawB=df[dexpression,2],
                       normA=x[dexpression,1],
                       normB=x[dexpression,2],
                       lg2A=log2(x[dexpression,1]),
                       lg2B=log2(x[dexpression,2]),
                       diff=(x[,2]-x[,1])[dexpression])
  
  
  ## make a copy of the results table ti display on the graph.
  display = results
  display$normA = format(display$normA, digits=2, nsmall=2)
  display$normB = format(display$normB, digits=2, nsmall=2)
  display$lg2A = format(display$lg2A, digits=2, nsmall=2)
  display$lg2B = format(display$lg2B, digits=2, nsmall=2)
  display$diff = format(display$diff, digits=2, nsmall=2)
  
  title=paste0(datasets,collapse=" vs. ")
  
  ## nf = layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(2,1), heights=c(1,1))
  
  ma.plot( rowMeans(log2(x)), log2(x[, 1])-log2(x[, 2]),
           xlab="Mean",ylab=paste0(datasets, collapse=" - "), cex=0.7) 
  textxy(rowMeans(log2(x))[dexpression],(log2(x[, 1])-log2(x[, 2]))[dexpression],genenames[dexpression])
  
  ## plot.new()
  ## addtable2plot(0,0,display,bty="o",display.rownames=TRUE,hlines=TRUE,
  ##               xpad=.1, ypad=.7, title="The table", cex=0.7)
  
#   ## draw the same plot to a PDF file.
#   postscript(file=paste0("/tmp/",title,".eps"), onefile=FALSE, horizontal=TRUE)
#   
#   ##    nf = layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(2,1), heights=c(1,1))
#   
#   ma.plot( rowMeans(log2(x)), log2(x[, 1])-log2(x[, 2]),
#            xlab="Mean",ylab=paste0(options, collapse=" - "), cex=0.7) 
#   textxy(rowMeans(log2(x))[dexpression],(log2(x[, 1])-log2(x[, 2]))[dexpression],genenames[dexpression])
#   
#   ##    plot.new()
#   ##    addtable2plot(0,0,display,bty="o",display.rownames=TRUE,hlines=TRUE,
#   ##                  xpad=.4, ypad=1, title="The table", cex=0.7)
#   dev.off()
  
  
  return (results)
}

