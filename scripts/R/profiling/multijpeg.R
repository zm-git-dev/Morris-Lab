

#gene="NR_029560"   # Mir150
gene='NM_001005419'	# 'Ado'  single exon reverse
gene='NM_013541'	# 'Gstp1'
gene='NM_011434'	# 'Sod1'
gene='NM_144903'	# Aldob
gene='NM_011044'	# 'Pck1'
gene="NM_007409"   # RPL22
gene='NM_016978'	# 'Oat'
gene="NM_007409"   # ADH1

genes = c('NM_001005419', 'NM_013541', 'NM_011434', 'NM_144903', 'NM_016978', "NM_007409")

for (gene in genes) {
    df = morris.getalignments("113010_A", gene)
        
    kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
    rownames(kg) <- kg$name
    gobj = transcript(kg[gene,])

    jpeg(filename = paste0("Rplot_", gobj$name2(), ".jpeg"), width=1200, height=909,type="quartz")
    
    par(mar=c(1,2,1,2))
    profile(df, gobj, bias="middle", minlen=28)
    dev.off()
}

