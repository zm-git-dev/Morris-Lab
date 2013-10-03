## plot multiple overlapping densvity functions for a single gene in
## multiple datasets.  this is one waypoint along with way to
## calculating density functions for small portions of the profile in
## order to find pause sites.  thrn we want to be able to compare
## pause site positions and profiles for a single gene across multiple
## experiments


require("sm")   ## for sm.density.compare


source("../morrislib.R")
source("profiling.R")

datasets= c("061113_A", "061113_B", "061113_C", "061113_D")
genome <- morris.getGenome(datasets[[1]])
descriptions = morris.fetchdesc(datasets)
motiflen = 9

profiles <- list(
    RPS23 = list(c()),
    RPS2 = list(c()),
    EEF2 = list(c()),
    EIF3E = list(c()),
    B2M = list(c())
)

for (gene in names(profiles)) {

    kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
    rownames(kg) <- kg$name
    stopifnot(nrow(kg) == 1)

    refseq = kg[1,'name']
    bdf = data.frame()

    for (dataset in datasets) {
        df = morris.getalignments(dataset, refseq)
        attr(df, "dataset") <- descriptions[dataset,"description"]
        ribo.profile = profile(df, kg[refseq,])
        print(paste0(refseq,":", ribo.profile$transcript()$name2()))
        df <- ribo.profile$plotpositions()

        df$exp = match(dataset,datasets)
        bdf <- rbind(bdf, df)

    }  ## for each dataset

    ## create value labels
    vlabels = factor(bdf$exp, levels=datasets, labels=descriptions$description)

    ## plot densities
    sm.density.compare(bdf$rpositions, bdf$exp, model="equal", xlab="position on gene")
    title(main=paste0(gene, " Ribosome density"))

    colfill<-c(2:(2+length(levels(vlabels)))) 
    legend(locator(1), levels(vlabels), fill=colfill)
    locator(1)


} ## for-each gene
