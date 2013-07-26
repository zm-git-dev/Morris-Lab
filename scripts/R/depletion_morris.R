library(IRanges)

datasets=c("030713_A.rrna.dist", "032513_A.rrna.dist")
dslist <- lapply(datasets, function(ds) read.table(paste0("~/Morris-Lab/analysis/030713/rrna_alignments/",ds), quote="\""))
stopifnot(all(unlist(lapply(dslist, function(ds) (all(c("V1","V2","V3","V4")== names(ds)))))))

df = merge(dslist[[1]],dslist[[2]],  by="V4",suffixes=c(".A",".B"))

## remove entries that have fewer than mincount reads.
mincount = 500
row.sub = apply(df[c("V1.A", "V1.B")], 1, function(row) (all(!is.na(row)) && all(row >= mincount)))
df = df[row.sub,]

## scaling the counts of rrna abundance
df[,c("V1.A", "V1.B")] = apply(df[,c("V1.A", "V1.B")], 2, function(col) col/sum(col))

IR_rna=IRanges(df$V3.A,df$V3.A+nchar(as.character(df$V4)))

## plot the points that would be seleted by morris oligos.
par(fig=c(0,0.50,0.1,0.9))
par(mgp=c(2,1,0))
plot(df$V1.A,df$V1.B,log="xy",cex=.6,xlab=datasets[1],ylab=datasets[2])

oligos <- read.delim("~/Morris-Lab/genome/mm9/oligos.bed", header=F, sep="",comment.char="#")
stopifnot(length(oligos) > 0)
stopifnot(all(c("V1","V2","V3","V4") == names(oligos)))
IR_oligo=IRanges(oligos$V2,oligos$V3)
countOverlaps(IR_rna, IR_oligo, minoverlap=20)
ov <- countOverlaps(IR_rna, IR_oligo, minoverlap=20)>0
points(df[ov,"V1.A"],df[ov,"V1.B"],pch=16,cex=1.1,col="red")
title("Depletion by Morris oligos", sub="scaled to by total fragments", cex.sub=.7)


par(fig=c(0.5,1,0.1,0.9), new=TRUE)
##par(mgp=c(1,1,0))
plot(df$V1.A,df$V1.B,log="xy",cex=.6,xlab=datasets[1],ylab=datasets[2])

oligos <- read.delim("~/Morris-Lab/genome/mm9/ingolia.bed", header=F, sep="",comment.char="#")
stopifnot(length(oligos) > 0)
stopifnot(all(c("V1","V2","V3","V4") == names(oligos)))
IR_oligo=IRanges(oligos$V2,oligos$V3)
ov <- countOverlaps(IR_rna, IR_oligo, minoverlap=5)>0
points(df[ov,"V1.A"],df[ov,"V1.B"],pch=16,cex=1.1,col="blue")
title("Depletion by Ingolia oligos", sub="scaled to by total fragments", cex.sub=.7)