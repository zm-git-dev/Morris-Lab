source("../morrislib.R")
source("profiling.R")

datasets= c("061113_A", "061113_B", "061113_C", "061113_D")

profiles <- list(
    RPS2 = list(c()),
    RPS23 = list(c()),
    EIF3E = list(c()),
    EIF3H = list(c()),
    EIF3I = list(c()),
    EEF1A1 = list(c()),
    EIF4EBP2 = list(c())
    
)

# prostate datasets
datasets <- c("103112_B_MM1",
              "030713_A_MM1",
              "032513_B_MM1",
              "103112_A_MM1",
              "030713_B_MM1",
              "032513_A_MM1")
profiles <- list(
    ADH1 = list(c())
)

genome <- morris.getGenome(datasets[1])
descriptions = morris.fetchdesc(datasets)


main <- function() {

    filename = "profiles.pdf"
    ##pdf(file=filename)
    par(mfrow = c(2,2))

    for (gene in names(profiles)) {
        kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
        rownames(kg) <- kg$name
        stopifnot(nrow(kg) == 1)
        
        refseq = kg[1,'name']

        for (dataset in datasets) {
            df = morris.getalignments(dataset, refseq)
            ##attr(df, "dataset") <- descriptions[dataset,"description"]
            ribo.profile = profile(df, kg[refseq,])
            print(paste0(refseq,":", ribo.profile$transcript()$name2()))
            for (x in profiles[[gene]]) {
                ## set the margin to give a more pleasing layout.
                ##     c(bottom, left, top, right) gives the number of
                ##     lines of margin to be specified on the four sides of the plot
                par(mar=c(1,2,1,2))
                print(x)
                plot(ribo.profile, minlen=28, xlim=x, units="aa")
            }  ## for each gene section
        }  ## for each dataset
    } ## for-each gene
    ##dev.off()
}

main()
