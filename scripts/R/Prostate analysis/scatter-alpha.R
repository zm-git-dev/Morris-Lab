

## try to reproduce a scatter plot like Dave showed.
## - eliminate all Riken genes
## - eliminate non-coding genes.
## - eliminate lowest expression TR- control sets. (<500K total reads)
## - eliminate the 3 least developed tumors from the TR+ sample set

source("../morrislib.R")

library(Biobase)
require(genefilter)


tumorSpecific <- c( "Tubb5", "Hist2h2bb", "Nov", "Chgb", "Hist1h1a" )
epithelial <- c("Tgm4", "Pbsn", "Sbp")
fibroblast <- c("Vim", "Itgb1", "Itga1", "Col1a1", "Col1a2")

# divide the data into two sets - control and treated...
##treated <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+/TR+")
treated <- c("103112_A", "030713_B", "041713_B")
control <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+")
datasets <- c(treated,control)
df.raw = morris.genecounts(datasets, group=NULL)
genome = attr(df.raw, "genome")

kg <- morris.getknowngenes(genome)

## eliminate lowest expression TR- control sets. (<600K total reads)
control <- control[colSums(df.raw[control], na.rm = TRUE) > 600000]

df <- df.raw[, c(control, treated)]

## Add a pseudocount to all missing data so logarithms taken later
## will be well-defined.
df[is.na(df)] <- 1

## toss out non-coding genes
df <- df[grep("NR_", rownames(df), invert=TRUE),]

## toss RIKEN genes.
df <- df[grep("Rik", kg$name2[match(rownames(df), kg$name)], invert=TRUE), ]

## toss any gene that does not have a corresponding common name
df <- df[!is.na(kg$name2[match(rownames(df), kg$name)]), ]

## scale all gene counts to the count of Col1a2
Col1a2 = with(kg, name[match("Col1a2", name2)])
n = as.numeric(df[Col1a2, control])
df[control] = sweep(df[control], 2, n, '/')
df[control] = sweep(df[control], 2, mean(n), '*')
n = as.numeric(df[Col1a2, treated])
df[treated] = sweep(df[treated], 2, n, '/')
df[treated] = sweep(df[treated], 2, mean(n), '*')

## eliminate any gene that has a count < 20
df <- df[apply(df, 1, function(row) (all(row > 20))),]

## Normalize to mapped reads per million mapped reads
attr(df, "genome") = genome
df.norm = morris.normalize(df[treated], normalization="rpm")
df.norm[control] = morris.normalize(df[control], normalization="rpm")

## Save the refseq names of fibroblast genes that we will want to highlight
rows.epithelial <- with(kg, name[match(epithelial, name2)])
rows.fibroblast <- with(kg, name[match(fibroblast, name2)])
rows.tumorSpecific <- with(kg, name[match(tumorSpecific, name2)])
rows.neither = rownames(df.norm)[(!rownames(df.norm) %in% c(rows.fibroblast, rows.epithelial, rows.tumorSpecific))]


m1=log2(rowMeans(df.norm[,control]))
m2=log2(rowMeans(df.norm[,treated]))

xlim=c(min(m1), max(m1))
ylim=c(min(m2), max(m2))
plot(m1[rows.neither], m2[rows.neither],
     xlab="Control", ylab="Tramp+",
     cex=.7, xlim=xlim, ylim=ylim, pch=20, col="mediumblue")


tstr = paste("Prostate Control (", length(control), ") vs TRAMP+ (", length(treated), ")")
ststr = "prostate Analysis/scatter.R, rpm,  norm(Col1a2), mean vs. mean"
title(tstr, sub=ststr, cex.main=1, cex.sub=1, col.sub="grey73")
    
points(m1[rows.epithelial], m2[rows.epithelial], pch=18, cex=1.2, col="red")		## diamond
points(m1[rows.fibroblast], m2[rows.fibroblast], pch=17, cex=1.2, col="yellow")		## triangle
points(m1[rows.tumorSpecific], m2[rows.tumorSpecific], pch='+', cex=1.2, col="green")	## square

# standard line of best fit - black line
abline(lm(m2[rows.fibroblast] ~ m1[rows.fibroblast]))

# force through [0,0] - blue line
##abline(lm(y ~ x + 0, data=test), col="blue")
abline(lm(m2 ~ m1), col="blue")

stopifnot(length(rows.epithelial) == length(epithelial))
stopifnot(length(rows.fibroblast) == length(fibroblast))
stopifnot(length(rows.tumorSpecific) == length(tumorSpecific))

repeat {
    ans <- identify(m1, m2, n=1, plot=FALSE)
    if(!length(ans)) break
    print(paste0("(",m1[ans], ",", m2[ans],")  ",rownames(df.norm)[ans]," ", with(kg, name2[match(rownames(df.norm)[ans], name)])))
}

Sys.sleep(0) 
s1 = rowSds(df.raw[,control])
s2 = rowSds(df.raw[,treated])
ttest=(m2-m1)/sqrt(s1^2/length(control)+s2^2/length(treated))
hist(ttest,nclass=100)

Sys.sleep(0) 
## Calculate the two-tailed T-test statistic
## probabiliy of seeing the a value greater than the measured value on 
## both the positive and negative extremese.
pval=2*(1-pt(abs(ttest),4))
hist(pval, nclass=100)
Sys.sleep(0) 
