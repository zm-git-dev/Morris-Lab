

#gene="NR_029560"   # Mir150
gene='NM_001005419'	# 'Ado'  single exon reverse
gene='NM_144903'	# Aldob
gene='NM_011044'	# 'Pck1'
gene="NM_007409"   # RPL22
gene='NM_016978'	# 'Oat'
gene='NM_011434'	# 'Sod1'
gene='NM_013541'	# 'Gstp1'
gene="NM_007409"   # ADH1

# read all the alignments from 'dataset' that align to 'gene'
# so these will be multiple alignments, but only those that align to
# a single gene.
#
plotIDs <- matrix(c(1:4), 2, 2, byrow=T)
layout(plotIDs, widths = rep(1.0/ncol(plotIDs), ncol(plotIDs)), heights =  rep(1.0/nrow(plotIDs), nrow(plotIDs)))
layout.show(4)

par(oma=c(1,1,1,1))

for (dataset in c("041713_A", "041713_B", "032513_B", "032513_A")) {
    
    df = morris.getalignments(dataset, gene)
    
    kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
    rownames(kg) <- kg$name
    gobj = transcript(kg[gene,])
    
    par(mar=c(1,2,1,2))
    profile(df, gobj, bias="middle", minlen=28)
}
