

df = morris.getalignments("113010_A")

kg <- morris.getknowngenes(attr(df, "genome"), group=NULL)
##kg <- data.frame(genome="mm9", name=c("gene1","gene2"), strand='+', txStart=1, txEnd=20,
##                 cdsStart=c(5,6), cdsEnd=c(8,9), exonCount=1, exonStarts="1,", exonEnds="20,",
##                 name2=c("test","test2"),stringsAsFactors = F)
rownames(kg) <- kg$name
##print(kg[gene,])

#df <- data.frame(position=c(49562310,49562310), length=22)
#attr(df,"genome") <- "mm9"

for (gene in rownames(df)) {
    print(profile(df, kg[gene,])$plotpositions())
}
