


profiles = list(
Alb=list(list(1,100)),
Fgg=list(list(1,100)),
Fgb=list(list(1,100)),
Apob=list(list(1,100)),
Apoa1=list(list(1,100))
)


for (gene in names(profiles)) {
    kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
    rownames(kg) <- kg$name
    stopifnot(nrow(kg) == 1)

    refseq = kg[1,'name']
    df = morris.getalignments("113010_A", refseq)

    ribo.profile = profile(df, kg[refseq,])
    print(paste0(refseq,":", ribo.profile$transcript()$name2()))
    for (x in profiles[[gene]]) {
        xlim=as.numeric(x)
        jpeg(filename = paste0(ribo.profile$transcript()$name2(), "_", xlim[1], "-", xlim[2], ".jpeg"),
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

