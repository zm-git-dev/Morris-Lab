#!/usr/bin/env Rscript

source("morrislib.R")
library(optparse) 

option_list <- list(
  make_option("--mincount", type="integer", default=25,
              help="Minimum number of reads for a gene [default %default]"),
  make_option("--normalization", default="quantile",
              help = "Function to normalize data, \"quantile\" or \"scale\" [default \"%default\"]"),
  make_option("--group", default=NULL,
              help="Database access group, usually \"morrisdata\" or \"remote\" [default %default]")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list),positional_arguments = TRUE)

if (length(opt$args) == 0) {
    ## If the user did not supply any arguments (running interactively
    ## inside RStudio), then prompt for the names of datasets.
    message("Type the names of datasets to compare, seperated by spaces")
    datasets = scan(what="list",nlines=1)
} else {
    ## otherwise the dataset names were supplied on the commandline.
    datasets = opt$args
}

print(opt$options$group)

results = morris.maplot(datasets, mincount=opt$options$mincount,
                       normalization=opt$options$normalization, group=opt$options$group)
print(results, digits=3)

message("Press Return To Continue")
invisible(readLines(n=1))
