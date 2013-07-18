#!/usr/bin/env Rscript

## source a file from wherevere this script was found.
## http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script/1815743
sourcelocal <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}

sourcelocal("morrislib.R")
library(optparse) 

option_list <- list(
  make_option("--mincount", type="integer", default=NA,
              help="Minimum number of reads for a gene [default %default]"),
  make_option("--normalization", default=NULL,
              help = "Function to normalize data, \"quantile\" or \"scale\" [default \"%default\"]"),
  make_option("--group", default=NULL,
              help="Database access group, usually \"morrisdata\" or \"remote\" [default %default]"),
  make_option("--report", type="integer", default=5,
              help="Report and label the N-most differentially expressed genes [default %default]"),
  make_option("--list", action="store_true", default=FALSE,
              help="Show the names of available datasets and exit.")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list),positional_arguments = TRUE)
if (opt$options$list) {
    print("Available datasets:")
    print(morris.datasets(group=opt$options$group));
} else {
    if (length(opt$args) == 0) {
        ## If the user did not supply any arguments (running interactively
        ## inside RStudio), then prompt for the names of datasets.
        message("Type the names of datasets to compare, seperated by spaces")
        datasets = scan(what="list",nlines=1)
    } else {
        ## otherwise the dataset names were supplied on the commandline.
        datasets = opt$args
    }

    group=opt$options$group
    normalization=opt$options$normalization
    mincount = opt$options$mincount
    
    df = morris.genecounts(datasets, group=group)

    ## Identify entries that do not have proper read depth, but don't
    ## eliminate them yet.  Some normalization methods rely upon
    ## having the full dataset available, even if only a portion of
    ## the dataset will eventually be reported.
    ##
    row.sub = NA
    if (!is.na(mincount)) 
      row.sub = apply(df, 1, function(row) (all(!is.na(row)) && all(row >= mincount)))

    if (!is.null(opt$options$normalization)) {
        x <- morris.normalize(df, normalization=normalization, group=group)
    

        ## merge the raw and normalized values into a single dataframe
        df <- merge(df, x, by=0,
                    suffixes=c("", paste0(".", normalization)), all=TRUE)
        rownames(df) <- df$Row.names
        df <- subset(df, select=-c(Row.names) )
    }

    ## Now remove the entries identified earlier that do not have
    ## proper read depth.
    ##
    if (!is.na(row.sub))
      df <- df[row.sub,, drop=FALSE]
    
    ## add common names for the refseq genes.
    df$common <- morris.commonnames(rownames(df), group=group)

    ## rearrange the columns to move the last column (common name) to the second column
    ## this is a little tricky b/c number of cols in df is variable.
    df <- df[c(ncol(df), 1:(ncol(df)-1))]

    options("scipen"=100, "digits"=4)
    print(df)
}
