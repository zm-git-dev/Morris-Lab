


profiles = list(
  'NM_009654', ## Alb
  'NM_009692', ## Apoa1 apolipoprotein A-I
  'NM_009693', ## Mus musculus apolipoprotein B (Apob)
  'NM_181849', ## FGB
  'NM_133862'  ## Fgg
  
)


for (gene in profiles) {
    df = morris.getalignments("113010_A", gene)
    kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
    rownames(kg) <- kg$name

    ribo.profile = profile(df, kg[gene,])
    print(paste0(gene," : ", ribo.profile$transcript()$name2() ))
    jpeg(filename = paste0("Rplot_", ribo.profile$transcript()$name2(), ".jpeg"),
         width=1200, height=909,type="quartz")

    ## set the margin to give a more pleasing layout.
    ##     c(bottom, left, top, right) gives the number of
    ##     lines of margin to be specified on the four sides of the plot
    par(mar=c(1,2,1,2))
    plot(ribo.profile, minlen=28, units="aa")
    dev.off()

}

