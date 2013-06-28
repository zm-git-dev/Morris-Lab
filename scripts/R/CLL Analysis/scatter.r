library(Biobase)
require(genefilter)

epithelial <- c("9530053A07Rik", "Tgm4", "Pbsn", "Rnase1", "Sbp")
fibroblast <- c("Vim", "Itgb1", "Itga1", "Col1a1", "Col1a2")

## treated <- morris.datasets(organism="Human", tissue="CLL", genotype="Ribo+/Col+")
treated <- c("061113_B")
control <- c("061113_A")
datasets <- c(treated, control)


ds = morris.genecounts(datasets, group=NULL)
genome = attr(ds, "genome")

df = ds
## Get rid of any rows containing NA.
df = df[complete.cases(df),]

## Exclude genes that are highly expressed in epithelial genes.
##cn <- morris.commonnames(rownames(df), genome, NULL) 
##df = df[!(cn %in% epithelial),]

## Exclude genes that have too few reads.
## row.sub = apply(df[,grep("103112",names(df),value=TRUE),drop = FALSE], 1, function(row) (any(row >= 30)))
## row.sub = row.sub & apply(df[,grep("030713",names(df),value=TRUE),drop = FALSE], 1, function(row) any(row >= 50))
## row.sub = row.sub & apply(df[,grep("032513",names(df),value=TRUE),drop = FALSE], 1, function(row) any(row >= 100))
## df = df[row.sub,]

## Normalize to mapped reads per million mapped reads
##  http://stackoverflow.com/questions/13830979/#13831155
df = sweep(df,2,(colSums(df)/(10^6)), '/')

## ## scale all RPM count to the RPM count of Col1a2
## ## 	First lookup the values for Col1a2
## cn <- morris.commonnames(rownames(df), genome, NULL) 
## reference = df[match("Col1a2", cn),]
## ##	Now divide each column in each row by to corresponding column of the reference.
## ##	http://stackoverflow.com/questions/13830979/#13831155
## df <- sweep(df,2,as.numeric(reference), '/')

## ## Save the refseq names of fibroblast genes that we will want to highlight
## hilite=names(cn[match(fibroblast,cn)])

e=(df)
m1=e[,control]
m2=e[,treated]
plot(m1, m2, xlab=paste0(control, collapse="-"), ylab=paste0(treated, collapse="-"),log="xy")
tstr = paste("control (", length(control), ") vs treated (", length(treated), ")")
title(tstr, sub="min=50, rpm",cex.main=.7,
      cex.sub=.6, col.sub="grey73")

