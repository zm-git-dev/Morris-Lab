#!/usr/bin/env Rscript

suppressMessages( library(RMySQL) )
suppressMessages( library(calibrate) )
suppressMessages( library(affy) )
library(preprocessCore)

suppressMessages( require(maptools) )
suppressMessages( library(plotrix) )
library(optparse) 

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

fetchData <- function(datasets, mincount, group) {
    message("Fetching table from database (this could take a while...)")
    con = dbConnect(dbDriver("MySQL"), group=group)

    qdata <- paste0(
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
    df.data= dbGetQuery(con, qdata)

    qdesc <- paste0("select d.name,e.description from experiments_tbl e join datasets_tbl d on d.expr_id=e.id where ",
                    paste0("d.name like '", datasets, "'", collapse=" or "))
    print(qdesc)
    df.desc <- dbGetQuery(con, qdesc)
    rownames(df.desc) = df.desc[,1]
    df.desc <- subset(df.desc, select = -c(name) )
    
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

scatter <- function(datasets, mincount=25, group="morrisdata", normalization="quantile") {
    datalist <- fetchData(datasets, mincount, group)
    df = datalist$df
    df.desc = datalist$descriptions

    ## identify each row with the gene name.
    rownames(df) = df[,1]
    genenames = df[,2]
    df <- subset(df, select = -c(name,geneSymbol) )

    ## normalize to total reads in each experiment.
    x <- switch(normalization,
                quantile = normalize.quantiles(data.matrix(df[,1:2])),
                TC = (df/colSums(df))*mean(colSums(df)),
                stop("Unrecognized 'normalization' value: ", normalization))
        
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

    ## make a copy of the results table to display on the graph. Number are formatted here with only
    ## a few decimal places.
    display = results
    display$normA = format(display$normA, digits=2, nsmall=2)
    display$normB = format(display$normB, digits=2, nsmall=2)
    display$lg2A = format(display$lg2A, digits=2, nsmall=2)
    display$lg2B = format(display$lg2B, digits=2, nsmall=2)
    display$diff = format(display$diff, digits=2, nsmall=2)
    
    title=paste0(datasets,collapse=" vs. ")

    ## nf = layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(2,1), heights=c(1,1))

    ##  Plot the data, hiding the points for now to prevent the call to grid()
    ##  and abline() from drawing over the points.
    options(scipen=100)
    plot(x[,1], x[,2], type = "n", log = "xy", main = "Log-log Plot",
         xlab=df.desc[datasets[1],1], ylab=df.desc[datasets[2],1], cex=0.7)
    grid(col = colors()[ 440 ], equilogs = FALSE)
    
    points( x[,1], x[,2])
    abline(lm(x[,2]~x[,1]), untf=TRUE, col="red")
    textxy(x[dexpression,1], x[dexpression,2], genenames[dexpression])

    ## plot.new()
    ## addtable2plot(0,0,display,bty="o",display.rownames=TRUE,hlines=TRUE,
    ##               xpad=.1, ypad=.7, title="The table", cex=0.7)
    
    ## draw the same plot to a PDF file.
    postscript(file=paste0("/tmp/",title,".eps"), onefile=FALSE, horizontal=TRUE)

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
results = scatter(opt$args, mincount=opt$options$mincount, normalization=opt$options$normalization, group=opt$options$group)
print(results, digits=3)

message("Press Return To Continue")
invisible(readLines("stdin", n=1))


