library(qarp)

## open a datafile containing read depth data
data <- qarp.read.cachedtsv("inst/extdata/PEO/PEO1-RPT-top200-coverage.tsv.gz")

## Establish which columns in the data frame contain read norm data
## We use columns of data in which each read depth has been normalized to the
## total number of aligned reads in the sample.
control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
datasets = c(control = control, treated = treated)

# eliminate extra columns that we don't need.
data = data[, c("symbol", "transPos", "exonNumber", datasets)]
data$symbol <- as.factor(data$symbol)

## create a list of factors that will be used to divide the datasets into
## control and treated groups.
aVec <- as.factor(sapply(datasets %in% control, function(x) if (x) {"control"} else {"treated"}))
names(aVec) <- datasets

## Draw the profiles of the LCN2 gene without normalizing the reads in any way
## This gives a good representation of raw read depth across the gene.
mat <- qarp.matrix(data, "LCN2", aVec)
qarp.plotProfile(mat, aVec)

## Draw a multidimensional sampling scatter plot derived from the
## read depth data.
mat <- qarp.matrix(data, "LGALS3BP", aVec)
qarp.plotMDS(mat, aVec)

## calculate and print a pvalue for just the ACTB gene
df <- qarp(data, aVec, genelist=c("ACTB"))
print(df)

## calculate pvalue for all genes in the dataset (genomewide)
## Order the list according to pvalue, draw  histogram of the
## distribution, and print the list
df = qarp(data, aVec)
df <- df[order(df$pvalues),,drop=FALSE]
hist(df$pvalues, n=50, main=attr(data, "filename"))
print(df)

