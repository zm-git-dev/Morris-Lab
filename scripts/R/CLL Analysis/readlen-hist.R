## There is a prominent cluster of reads towards the 3'-end of RPS23.
## This cluster may be due to a ribosome pause site, or to a protein
## complex that co-pecipitates with other fragments inside of
## ribosomes.
##
## a working hyposesis is that the distribution of reads protected by
## ribosomes will be different than the distribution of fragments
## protected by a protein complex.


## get reads from RPS23.
## Consider only those reads that are within a certain range of the message.
## draw a histogram of their length distribution.

source("../morrislib.R")
source("../profiling/profiling.R")


datasets <- c("061113_A", "061113_B", "061113_C", "061113_D")
genome <- morris.getGenome(dataset[1])
descriptions <- morris.fetchdesc(datasets)

gene <- "RPS23"

kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
rownames(kg) <- kg$name
stopifnot(nrow(kg) == 1)

refseq <- kg[1,'name']

clusterA <- c()
clusterB <- c()
for (dataset in datasets) {
    df = morris.getalignments(dataset, refseq)
    attr(df, "dataset") <- descriptions[dataset,"description"]
    ribo.profile = profile(df, kg[refseq,])
    print(paste0(refseq,":", ribo.profile$transcript()$name2()))
    x=NULL

    ## set the margin to give a more pleasing layout.
    ##     c(bottom, left, top, right) gives the number of
    ##     lines of margin to be specified on the four sides of the plot
    par(mar=c(1,2,1,2))
    df <- plot(ribo.profile, minlen=28, xlim=x, units="aa")

    x <- with(df, df[rposition>=95 && rposition<=100, 'length'])
    clusterA = c(clusterA, x)

    x <- with(df, df[rposition < 95,'length'])
    clusterB = c(clusterB, x)
}


