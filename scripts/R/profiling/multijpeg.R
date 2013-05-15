

#gene="NR_029560"   # Mir150
gene='NM_001005419'	# 'Ado'  single exon reverse
gene='NM_013541'	# 'Gstp1'
gene='NM_011434'	# 'Sod1'
gene='NM_144903'	# Aldob
gene='NM_011044'	# 'Pck1'
gene="NM_007409"   # RPL22
gene='NM_016978'	# 'Oat'
gene="NM_007409"   # ADH1

profiles = list(
NM_007409=list(list(15,115), list(90,190), list(176,276), list(220,320), list(260,360)), 
NM_011434=list(list(15,115)),
NM_013541=list(list(1,80)),
NM_011044=list(list(30,130), list(100,200), list(160,260), list(520,620))
)


for (gene in names(profiles)) {
    df = morris.getalignments("113010_A", gene)
    kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
    rownames(kg) <- kg$name

    ribo.profile = profile(df, kg[gene,])
    print(paste0(gene,":"))
    for (x in profiles[[gene]]) {
        xlim=as.numeric(x)
        jpeg(filename = paste0("Rplot_", gobj$name2(), "_", xlim[1], "-", xlim[2], ".jpeg"),
             width=1200, height=909,type="quartz")

        ## set the margin to give a more pleasing layout.
        ##     c(bottom, left, top, right) gives the number of
        ##     lines of margin to be specified on the four sides of the plot
        par(mar=c(1,2,1,2))
        print(xlim)
        plot(ribo.profile, minlen=28, xlim=xlim, units="aa")
        dev.off()
    }

}

