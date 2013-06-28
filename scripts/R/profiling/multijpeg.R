

dataset="113010_A"
genome = morris.getGenome(dataset)

profiles = list(
  Actb=list(c(), c(1,250)),
  Cat=list(c(), c(1,250), c(225,350), c(325)),
  Gnmt=list(c(), c(1,150))
)

main <- function() {

    for (gene in names(profiles)) {
        kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
        rownames(kg) <- kg$name
        stopifnot(nrow(kg) == 1)
        
        refseq = kg[1,'name']
        df = morris.getalignments(dataset, refseq)
        
        ribo.profile = profile(df, kg[refseq,])
        print(paste0(refseq,":", ribo.profile$transcript()$name2()))
        for (x in profiles[[gene]]) {
            xlim=as.numeric(x)

            ## Create a filename to store the jpeg image in
            suffix = ""
            if (length(xlim) == 2)
              suffix=paste0( "_", xlim[1], "-", xlim[2])
            else if (length(xlim) == 1) {
                suffix=paste0( "_", xlim[1], "-")
            }
            
            filename = paste0(ribo.profile$transcript()$name2(), suffix, ".jpeg")

            jpeg(filename=filename, width=1200, height=909,type="quartz")

            ## set the margin to give a more pleasing layout.
            ##     c(bottom, left, top, right) gives the number of
            ##     lines of margin to be specified on the four sides of the plot
            par(mar=c(1,2,1,2))
            print(x)
            plot(ribo.profile, minlen=28, xlim=x, units="aa")
            dev.off()
        }
    }
}

main()
