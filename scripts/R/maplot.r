#!/usr/bin/env Rscript

suppressMessages( library(RMySQL) )
suppressMessages( library(calibrate) )
suppressMessages( library(affy) )
library(preprocessCore)

suppressMessages( require(maptools) )
suppressMessages( library(plotrix) )
library(optparse) 

buildQuery <- function(datasets, mincount) {
  query <- paste0(
    'select name, geneSymbol,\n',
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
  ##	, paste0(" limit 100; "))
  
  return (query)
}

## fetchData - retrieve a dataframe containing the data we are interested in comparing.
## in this case we are retrieving the data from a sql database.
##
## options - contains the names of two datasets, e.g. c("110112_A_MM1", "110112_B_MM1")
##
## mincount - the minimum number of reads that a gene expression must possess before it will be counted.
##
## group - The mysql database group to use.  "morrisdata" assumes the
## database is accessible on the local network.  "remote" assumes the
## database is accessible on localhost:6607.  Usually this means you
## are at a remote location and you have set up a tunnel via ssh to
## the database at UW.

fetchData <- function(options, mincount, group) {
    message("Fetching table from database (this could take a while...)")
    query = buildQuery(options, mincount=mincount)
    con = dbConnect(dbDriver("MySQL"), group=group)
    df = dbGetQuery(con, query)

}

##
## options - contains the names of two datasets, e.g. c("110112_A_MM1", "110112_B_MM1")
##
## mincount - the minimum number of reads that a gene expression must possess before it will be counted.
##
## group - The mysql database group to use.  "morrisdata" assumes the
## database is accessible on the local network.  "remote" assumes the
## database is accessible on localhost:6607.  Usually this means you
## are at a remote location and you have set up a tunnel via ssh to
## the database at UW.

maplot <- function(options, mincount=25, group="morrisdata", normalization="quantile") {
    df <- fetchData(options, mincount, group)

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

    quartz()

    ## make a copy of the results table ti display on the graph.
    display = results
    display$normA = format(display$normA, digits=2, nsmall=2)
    display$normB = format(display$normB, digits=2, nsmall=2)
    display$lg2A = format(display$lg2A, digits=2, nsmall=2)
    display$lg2B = format(display$lg2B, digits=2, nsmall=2)
    display$diff = format(display$diff, digits=2, nsmall=2)
    
    title=paste0(options,collapse=" vs. ")

    ## nf = layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(2,1), heights=c(1,1))

    ma.plot( rowMeans(log2(x)), log2(x[, 1])-log2(x[, 2]),
            xlab="Mean",ylab=paste0(options, collapse=" - "), cex=0.7) 
    textxy(rowMeans(log2(x))[dexpression],(log2(x[, 1])-log2(x[, 2]))[dexpression],genenames[dexpression])

    ## plot.new()
    ## addtable2plot(0,0,display,bty="o",display.rownames=TRUE,hlines=TRUE,
    ##               xpad=.1, ypad=.7, title="The table", cex=0.7)
    
    ## draw the same plot to a PDF file.
    postscript(file=paste0("/tmp/",title,".eps"), onefile=FALSE, horizontal=TRUE)

##    nf = layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(2,1), heights=c(1,1))
    
    ma.plot( rowMeans(log2(x)), log2(x[, 1])-log2(x[, 2]),
            xlab="Mean",ylab=paste0(options, collapse=" - "), cex=0.7) 
    textxy(rowMeans(log2(x))[dexpression],(log2(x[, 1])-log2(x[, 2]))[dexpression],genenames[dexpression])

##    plot.new()
##    addtable2plot(0,0,display,bty="o",display.rownames=TRUE,hlines=TRUE,
##                  xpad=.4, ypad=1, title="The table", cex=0.7)
    dev.off()

      
    return (results)
}


option_list <- list(
                    make_option("--mincount", type="integer", default=25,
                                help="Minimum number of reads for a gene [default %default]"),
                    make_option("--normalization", default="quantile",
                                help = "Function to normalize data, \"quantile\" or \"scale\" [default \"%default\"]"),
                    make_option("--group", default="morrisdata",
                                help="Database access group, usually \"morrisdata\" or \"remote\" [default %default]")
                    )

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list),positional_arguments = TRUE)
results = maplot(opt$args, mincount=opt$options$mincount, normalization=opt$options$normalization, group=opt$options$group)
print(results, digits=3)

message("Press Return To Continue")
invisible(readLines("stdin", n=1))


