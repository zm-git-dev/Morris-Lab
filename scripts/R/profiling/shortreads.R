plotIDs <- matrix(c(1:4), 2, 2, byrow=T)
layout(plotIDs, widths = rep(1.0/ncol(plotIDs), ncol(plotIDs)), heights =  rep(1.0/nrow(plotIDs), nrow(plotIDs)))
layout.show(4)

par(oma=c(1,1,1,1))

df = morris.getalignments(dataset, gene)

kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
rownames(kg) <- kg$name
gobj = transcript(kg[gene,])
    
par(mar=c(1,2,1,2))
profile(df, gobj, bias="middle", minlen=28)
