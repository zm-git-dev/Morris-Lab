library(Biobase)
require(genefilter)

epithelial <- c("9530053A07Rik", "Tgm4", "Pbsn", "Rnase1", "Sbp")
fibroblast <- c("Vim", "Itgb1", "Itga1", "Col1a1", "Col1a2")

# divide the data into two sets - control and treated...
treated <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+/TR+")
control <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+")
datasets <- c(treated,control)
df = morris.genecounts(datasets, group=NULL)
genome = attr(df, "genome")

## Get rid of any rows that are missing data...
## Hmmmmm.   Maybe should keep those and just set them to zero?
##
df[is.na(df)] <- 0
##df = df[complete.cases(df),]

## Identify genes that have too few reads.
row.sub = apply(df[,grep("103112",names(df),value=TRUE),drop = FALSE], 1, function(row) (any(row >= 30)))
row.sub = row.sub & apply(df[,grep("030713",names(df),value=TRUE),drop = FALSE], 1, function(row) any(row >= 50))
row.sub = row.sub & apply(df[,grep("032513",names(df),value=TRUE),drop = FALSE], 1, function(row) any(row >= 100))

cn <- morris.commonnames(rownames(df), genome, NULL) 
stopifnot(nrow(df)==nrow(cn))
stopifnot(rownames(df)==rownames(cn))

## scale all RPM count to the RPM count of Col1a2
## 	First lookup the values for Col1a2
reference = df[match("Col1a2", cn$common),]
stopifnot(!any(is.na(reference)))

##	Now divide each column in each row by to corresponding column of the reference.
##	http://stackoverflow.com/questions/13830979/#13831155
df <- sweep(df,2,as.numeric(reference), '/')

## Normalize to mapped reads per million mapped reads
attr(df,"genome") = genome
df = morris.normalize(df, normalization="rpm")

## Exclude genes that have too few reads.
df = df[row.sub,,drop=FALSE]
cn = cn[row.sub,,drop=FALSE]

## Save the refseq names of fibroblast genes that we will want to highlight
rows.epithelial=rownames(cn)[cn$common %in% epithelial]
rows.fibroblast=rownames(cn)[cn$common %in% fibroblast]
rows.neither=rownames(cn)[!cn$common %in% c(fibroblast, epithelial)]
    
m1=rowMeans(df[,control])
m2=rowMeans(df[,treated])

xlim=c(min(m1), max(m1))
ylim=c(min(m2), max(m2))
plot(m1, m2, xlab="Control", ylab="Tramp+",log="xy",cex=.8)
plot(m1[rows.neither], m2[rows.neither], xlab="Control", ylab="Tramp+",log="xy",cex=.7,
     xlim=xlim, ylim=ylim)
tstr = paste("control (", length(control), ") vs TRAMP+ (", length(treated), ")")
title(tstr, sub="min=50, rpm",cex.main=.7, cex.sub=1, col.sub="grey73")
    
points(m1[rows.epithelial], m2[rows.epithelial], pch=18, cex=1, col="red")
points(m1[rows.fibroblast], m2[rows.fibroblast], pch=17, cex=1, col="blue")


repeat {
    ans <- identify(m1, m2, n=1, plot=FALSE)
    if(!length(ans)) break
    print(paste0("(",m1[ans], ",", m2[ans],")  ",rownames(cn)[ans]," ", cn[ans,]))
}



s1 = rowSds(e[,control])
s2 = rowSds(e[,treated])
ttest=(m2-m1)/sqrt(s1^2/length(control)+s2^2/length(treated))
hist(ttest,nclass=100)

## Calculate the two-tailed T-test statistic
## probabiliy of seeing the a value greater than the measured value on 
## both the positive and negative extremese.
pval=2*(1-pt(abs(ttest),4))
hist(pval, nclass=100)
