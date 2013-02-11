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

    ## Build an SQL query string for retrieving alignment data from the database.
    ##
    dataQuery <- paste0('select name, geneSymbol,\n',
                    paste0("\tt", seq(datasets), ".rawcount",collapse=",\n"),
                    "\nfrom knowngenes_tbl\n",
                    paste0("left outer join ",
                           "(select FPKM, rawcount, gene_id \n",
                           " from expression_tbl join datasets_tbl ON expression_tbl.dataset_id = datasets_tbl.id \n",
                           " where datasets_tbl.name like '", datasets, "' and rawcount>", mincount,
                           " ) t",seq(datasets)," on t",seq(datasets),".gene_id = id\n",collapse="\n"),
                    paste0(" where ",
                           paste0("t",seq(datasets),".gene_id is not null",collapse=" and ")
                           ))

    ## Build an SQL query sring for retrieving the descriptions of
    ## experiments from the database.  These are usually used to add
    ## useful labels to axes on plots.
    descQuery <- paste0("select d.name,e.description from experiments_tbl e join datasets_tbl d on d.expr_id=e.id where ",
                    paste0("d.name like '", datasets, "'", collapse=" or "))

    result = tryCatch({
        drv = dbDriver("MySQL")
        if (is.null(group)) {
            message("group is NULL, using default parameters to connect to database...")
            con = dbConnect(drv, user="readonly", password="readonly", dbname="morris", host="morrislab.bchem.washington.edu")
        } else  {
            message("group is  NOT NULL...")
            con = dbConnect(drv, group=group)
        }
        df.data= dbGetQuery(con, dataQuery)

        df.desc <- dbGetQuery(con, descQuery)
        rownames(df.desc) = df.desc[,1]
        df.desc <- subset(df.desc, select = -c(name) )
    }, warning = function(w) {
        print(w$message)
    }, error = function(e) {
        print(e$message)
        print(e)
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        if (exists("drv")) 
          dbUnloadDriver(drv)
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
  df <- subset(df, select = -c(name,geneSymbol) )
  
  ## normalize to total reads in each experiment.
  if (normalization=="quantile") {
      x <- normalize.quantiles(data.matrix(df[,1:2]))
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

