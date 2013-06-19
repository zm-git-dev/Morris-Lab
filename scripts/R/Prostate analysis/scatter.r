library(Biobase)
require(genefilter)

datasets <- c("103112_B_MM1",
              "030713_A_MM1",
              "032513_B_MM1",
              "103112_A_MM1",
              "030713_B_MM1",
              "032513_A_MM1")

ds = morris.genecounts(datasets, group=NULL)
genome = attr(ds, "genome")

df = ds
## Get rid of any rows containing NA.
df = df[complete.cases(df),]

## Exclude genes that are highly expressed in epithelial genes.
exclude = c("9530053A07Rik", "Tgm4", "Pbsn", "Rnase1", "Sbp")
cn <- morris.commonnames(rownames(df), genome, NULL) 
df = df[!(cn %in% exclude),]

## Exclude genes that have too few reads.
row.sub = apply(df[,grep("103112",names(df))], 1, function(row) (any(row >= 30)))
row.sub = row.sub & apply(df[,grep("030713",names(df))], 1, function(row) any(row >= 50))
row.sub = row.sub & apply(df[,grep("032513",names(df))], 1, function(row) any(row >= 100))
df = df[row.sub,]

## Normalize to mapped reads per million mapped reads
##  http://stackoverflow.com/questions/13830979/#13831155
df = sweep(df,2,(colSums(df)/(10^6)), '/')

## scale all RPM count to the RPM count of Col1a2
## 	First lookup the values for Col1a2
cn <- morris.commonnames(rownames(df), genome, NULL) 
reference = df[match("Col1a2", cn),]
##	Now divide each column in each row by to corresponding column of the reference.
##	http://stackoverflow.com/questions/13830979/#13831155
df <- sweep(df,2,as.numeric(reference), '/')

e=(df)
m1=rowMeans(e[,1:3])
m2=rowMeans(e[,4:6])
plot(m1, m2, xlab="Control", ylab="Tramp+",log="xy")
title("Control vs. TRAMP+ (3 samples of each)", sub="min=50,epithel excluded, rpm",cex.main=.7,
      cex.sub=.6, col.sub="grey73")

s1 = rowSds(e[,1:3])
s2 = rowSds(e[,4:6])
ttest=(m2-m1)/sqrt(s1^2/3+s2^2/3)
hist(ttest,nclass=100)

## Calculate the two-tailed T-test statistic
## probabiliy of seeing the a value greater than the measured value on 
## both the positive and negative extremese.
pval=2*(1-pt(abs(ttest),4))
hist(pval, nclass=100)
