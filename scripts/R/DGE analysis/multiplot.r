resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}


datasets <- c("103112_B_MM1",
              "030713_A_MM1",
              "032513_B_MM1",
              "103112_A_MM1",
              "030713_B_MM1",
              "032513_A_MM1")
epithelial <- c("9530053A07Rik", "Tgm4", "Pbsn", "Rnase1", "Sbp")
fibroblast <- c("Vim", "Itgb1", "Itga1", "Col1a1", "Col1a2")


ds = morris.genecounts(datasets, group=NULL)

df = ds
## Get rid of any rows containing NA.
df = df[complete.cases(df),]

## Remove epithelial genes
exclude = c("9530053A07Rik", "Tgm4", "Pbsn", "Rnase1", "Sbp")
cn <- morris.commonnames(rownames(df), NULL) 
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
cn <- morris.commonnames(rownames(df), NULL)
names(cn)<-rownames(df)
reference = df[match("Col1a2", cn),]
##	Now divide each column in each row by to corresponding column of the reference.
##	http://stackoverflow.com/questions/13830979/#13831155
df <- sweep(df,2,as.numeric(reference), '/')

## Save the refseq names of fibroblast genes that we will want to highlight
hilite=names(cn[match(fibroblast,cn)])

par(mgp=c(1,0,0))
par(fig=c(0.30,0.70,0,0.5))
e=(df)
m1=rowMeans(e[,1:3])
m2=rowMeans(e[,4:6])
plot(m1, m2, xlab="Control", ylab="Tramp+", log="xy", cex=.4)
title("Control vs. TRAMP+\nMean", sub="min=50,epithel excluded, rpm",cex.main=2,
      cex.sub=.6, col.sub="grey43")
points(tmp[hilite,], pch=20, cex=1.4, col="red", log="xy")

par(fig=c(0,0.35,0.5,1), new=TRUE)
tmp <- df[,grep("103112",names(df))]
plot(tmp, cex=.4, xlab="Control", ylab="Tramp+",log="xy")
title("103112", cex.main=1)
points(tmp[hilite,], pch=20, cex=1.4, col="red", log="xy")

par(fig=c(0.32,0.67,0.5,1), new=TRUE)
tmp <- df[,grep("030713",names(df))]
plot(tmp, cex=.4, xlab="Control", ylab="Tramp+",log="xy")
title("030713", cex.main=1)
points(tmp[hilite,], pch=20, cex=1.4, col="red", log="xy")


par(fig=c(0.65,.98,0.5,1), new=TRUE)
tmp <- df[,grep("032513",names(df))]
plot(tmp, cex=.4, xlab="Control", ylab="Tramp+",log="xy")
title("032513", cex.main=1)
points(tmp[hilite,], pch=20, cex=1.4, col="red", log="xy")


