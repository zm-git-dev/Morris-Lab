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
df[is.na(df)] <- 0

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
df.norm <- morris.normalize(df, norm)

## Now remove the entries identified earlier that do not have
## proper read depth.
##
if (!is.na(row.sub))
    df.norm <- df.norm[row.sub,, drop=FALSE]

## Save the refseq names of  genes that we will want to highlight
kg = morris.getknowngenes(genome)
rownames(kg) <- kg$name
snorna = kg$name[grep("SNOR", kg$name2)]

m1 = log2(apply(df.norm[,control,drop=FALSE], 1, mean))
m2 = log2(apply(df.norm[,treated,drop=FALSE], 1, mean))

xlim=c(min(m1), max(m1))
ylim=c(min(m2), max(m2))

plot(m1[!(rownames(df.norm) %in% snorna)], m2[!(rownames(df.norm) %in% snorna)],
     xlab=paste0(control, collapse="-"), ylab=paste0(treated, collapse="-"),
     xlim=xlim, ylim=ylim, pch=16, cex=0.7)

points(m1[rownames(df.norm) %in% snorna], m2[rownames(df.norm) %in% snorna], col='blue', pch=17, cex=1)

repeat {
    ans <- identify(m1, m2, n=1, plot=FALSE)
    if(!length(ans)) break
    print(paste0("(",m1[ans], ",", m2[ans],")  ",rownames(df.norm)[ans]," ", kg[rownames(df.norm)[ans],"name2"]))
}


descriptions = morris.fetchdesc(datasets, group=group)
tstr = paste0(descriptions[,1], collapse="\n vs. ")


# tstr = paste("control (", length(control), ") vs treated (", length(treated), ")")
title(tstr, sub=paste0("SNORNA", norm, collapse=", "),cex.main=.9, cex.sub=.6, col.sub="grey73")
