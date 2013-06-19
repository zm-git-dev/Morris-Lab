library(IRanges)
rna=data.frame(
    V1=c(9),
    V2=c(50),
    V3="rna1")

oligo=data.frame(
    V1=c(3),
    V2=c(9),
    V3="oligo1")



IR_rna=IRanges(rna$V1, rna$V2)
IR_oligo=IRanges(oligo$V1, oligo$V2)
countOverlaps(IR_rna, IR_oligo,minoverlap=1)