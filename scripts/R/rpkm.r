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

fetchKnownGenes <- function(group) {
    message("Fetching table from database (this could take a while...)")
    con = dbConnect(dbDriver("MySQL"), group=group)
    qdata <- paste0("SELECT exonStarts,CONVERT(exonStarts USING utf8) FROM mm9.knownGene",
                      " limit 10")
    
    df.data= dbGetQuery(con, qdata)
    
    return ( df.data )
}



