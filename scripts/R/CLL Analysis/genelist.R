source("../morrislib.R")

options("max.print"=30)

group=NULL
treated <- c("061113_B")
control <- c("061113_D")
datasets <- c(control)

datasets <- c('061113_A', '061113_B','061113_C','061113_D')

df <- morris.genecounts(datasets, group=group)
genome <- attr(df, "genome")
stopifnot(!is.null(attr(df, "genome")))

## replace NA with 0
df[is.na(df)] <- 0
stopifnot(!is.null(attr(df, "genome")))

## find common gene names for all the refeseq names
df$common <- morris.commonnames(rownames(df), genome, group=group)
stopifnot(!is.null(attr(df, "genome")))

## identify which genes are noncoding 
kg = morris.getknowngenes(genome)
kg$coding = (kg$cdsStart != kg$cdsEnd)
df$coding = kg[match(rownames(df), kg$name), 'coding']
stopifnot(!is.null(attr(df, "genome")))

## rearrange the columns to move the last column (common name) to the second column
## this is a little tricky b/c number of cols in df is variable.

## Copy the attribute for the original data.frame if necessary
## see http://stackoverflow.com/a/10420036/1135316
## we do this to preserve the 'genome' attribute and any others that may be used in
## the future.   However we don't want to save the names of columns because those have
## just been rearranged.   It's complicated.

saveattr = attributes(df)
df <- df[c(ncol(df)-1,ncol(df),1:(ncol(df)-2))]
n = names(df)
mostattributes(df) = saveattr
names(df) = n
stopifnot(!is.null(attr(df, "genome")))

rpkm <- morris.normalize(df, normalization="rpkm", group=group)
rpkm$common = df[match(rownames(df),rownames(rpkm)),'common']
rpkm$coding = df[match(rownames(df),rownames(rpkm)),'coding']
rpkm <- rpkm[c(ncol(rpkm)-1,ncol(rpkm),1:(ncol(rpkm)-2))]

stopifnot(!is.null(attr(df, "genome")))

write.csv(rpkm, file = "CLL.RPKM.csv")
write.csv(df, file = "CLL.RAW.csv")
