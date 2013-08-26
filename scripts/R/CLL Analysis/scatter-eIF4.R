source("../morrislib.R")
options("max.print"=30)

library(Biobase)
require(genefilter)

norm = "rpm"
mincount = 30
highlight <- c("RPS2", "RPS23", "EIF3E", "EIF3H", "EIF3I", "EEF1a1", "EIF4EBP2")

## treated <- morris.datasets(organism="Human", tissue="CLL", genotype="Ribo+/Col+")
control <- c("061113_A")
treated <- c("061113_D")
datasets <- c(treated, control)


df <- morris.genecounts(datasets, group=NULL)
genome <- attr(df, "genome")

## Get rid of any rows containing NA.
df[is.na(df)] <- 0


## Identify genes that have too few reads.
row.sub = NULL
if (!is.na(mincount)) {
    ## row.sub is the list of rows to KEEP!
    row.sub = apply(df, 1, function(row) (all(!is.na(row)) && all(row >= mincount)))
    df <- df[row.sub,]
}

# toss out non-coding genes
df <- df[grep("NR_", rownames(df), invert=TRUE),]

## Normalize to RPM
df.norm <- morris.normalize(df, norm)

## Save the refseq names of  genes that we will want to highlight
cn <- morris.commonnames(rownames(df), genome, NULL) 
stopifnot(rownames(df)==rownames(cn))

## Identify entries that do not have proper read depth, but don't
## eliminate them yet.  Some normalization methods rely upon
## having the full dataset available, even if only a portion of
## the dataset will eventually be reported.
##
rows.highlight <- rownames(df.norm)[cn$common %in% highlight]
rows.nohighlight <- rownames(df.norm)[!(cn$common %in% highlight)]

m1 = log2(apply(df.norm[,control,drop=FALSE], 1, mean))
m2 = log2(apply(df.norm[,treated,drop=FALSE], 1, mean))

xlim=c(min(m1), max(m1))
ylim=c(min(m2), max(m2))

descriptions = morris.fetchdesc(datasets)


plot(m1[rows.nohighlight], m2[rows.nohighlight],
     xlab=descriptions[control,1], ylab=descriptions[treated,1],
     xlim=xlim, ylim=ylim, pch=16, cex=0.7)

points(m1[rows.highlight], m2[rows.highlight], col='red', pch=17, cex=1)
textxy(m1[rows.highlight], m2[rows.highlight], cn[rows.highlight,"common"], dcol="red")

tstr = paste0(descriptions[,1], collapse="\n vs. ")
title(tstr, sub=paste0("eIF4 highlight", norm, collapse=", "), cex.main=.9, cex.sub=.6, col.sub="grey73")


##legend(locator(), c("Protein Coding", "eIF4 pathway"), col = c("black", "red"), pch = c(16, 17), merge=TRUE, bg="gray95", text.font=2)


repeat {
    ans <- identify(m1, m2, n=1, plot=FALSE)
    if(!length(ans)) break
    print(paste0("(",m1[ans], ",", m2[ans],")  ",rownames(df.norm)[ans]," ", cn[ans,"common"]))
}


