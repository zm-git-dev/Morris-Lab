source("../morrislib.R")
options("max.print"=30)

library(Biobase)
require(genefilter)

norm = "rpm"
mincount = 30

## treated <- morris.datasets(organism="Human", tissue="CLL", genotype="Ribo+/Col+")
control <- c("061113_A")
treated <- c("061113_D")
datasets <- c(treated, control)


df <- morris.genecounts(datasets, group=NULL)
genome <- attr(df, "genome")

## Get rid of any rows containing NA.
#df <- df[complete.cases(df),]

## Identify entries that do not have proper read depth, but don't
## eliminate them yet.  Some normalization methods rely upon
## having the full dataset available, even if only a portion of
## the dataset will eventually be reported.
##
row.sub = NA
if (!is.na(mincount)) {
    ## row.sub is the list of rows to KEEP!
    row.sub = apply(df, 1, function(row) (all(!is.na(row)) && all(row >= mincount)))
}

## Normalize to mapped reads per million mapped reads
df <- morris.normalize(df, norm)

## Now remove the entries identified earlier that do not have
## proper read depth.
##
if (!is.na(row.sub))
    df <- df[row.sub,, drop=FALSE]

## ## Save the refseq names of fibroblast genes that we will want to highlight
## hilite=names(cn[match(fibroblast,cn)])
kg = morris.getknowngenes(genome)
sno = grep("SNOR", kg$name2)
notsno = grep("SNOR", kg$name2, invert=TRUE)

m = df[match(kg$name[notsno], rownames(df)),]
m1=m[, control]
m2=m[, treated]
plot(m1, m2, xlab=paste0(control, collapse="-"),
     ylab=paste0(treated, collapse="-"),log="xy", pch=16, cex=0.7)

m = df[match(kg$name[sno], rownames(df)),]
m1=m[, control]
m2=m[, treated]
points(m1, m2,col='blue', pch=17, cex=1)

descriptions = morris.fetchdesc(datasets, group=group)
tstr = paste0(descriptions[,1], collapse="\n vs. ")


# tstr = paste("control (", length(control), ") vs treated (", length(treated), ")")
title(tstr, sub=paste0("SNORNA", norm, collapse=", "),cex.main=.9, cex.sub=.6, col.sub="grey73")
