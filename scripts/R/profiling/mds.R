

# Experiment with Milti-Dimensional Scaling (MDS) to compare
# distribution of ribosomes on a single gene in multiple experiments.
#
# We are going to take the position data for a gene in an experiment
# and turn it into a row of a matrix.  We will use dist on that matrix
# to compute the complete distance matrix betweem expeirments.  Then
# we will use cmdscale() to compute a low-dimension field that can be
# used to measure how close or far from one another the ribosome
# profiles are.
#
#

source("../morrislib.R")
source("profiling.R")
require("ggplot2")

# CLL datasets
datasets= c("061113_A", "061113_B", "061113_C", "061113_D")
profiles <- list(
    RPS23 = list(c())
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


genome <- morris.getGenome(datasets[[1]])
descriptions = morris.fetchinfo(datasets)[,"description", drop=FALSE]
print(descriptions)
stats = morris.fetchstats(datasets)

main <- function() {

    for (gene in names(profiles)) {
        kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
        rownames(kg) <- kg$name
        stopifnot(nrow(kg) == 1)
        
        refseq = kg[1,'name']
        mat <- matrix(0, 0, 50)
        for (dataset in datasets) {
            df = morris.getalignments(dataset, refseq)
            attr(df, "dataset") <- descriptions[dataset,"description"]
            ribo.profile = profile(df, kg[refseq,])
            print(paste0(refseq,":", ribo.profile$transcript()$name2()))
            df <- ribo.profile$plotpositions()
            for (x in profiles[[gene]]) {
                if (is.null(x)) {
                    xlim <- c(1, ribo.profile$transcript()$txLength())
                } else {
                    xlim = as.numeric(x)
                }

                ## count how many reads occur on each position.
                histdata <- hist(df$rposition, breaks=c(1:ribo.profile$transcript()$txLength()), plot=FALSE)
                scores <- histdata$counts * (mean(stats[datasets, "raw_count"])/stats[dataset,"raw_count"])
                scores <- histdata$counts
                dim(mat) = c(dim(mat)[1], length(scores))
                mat <- rbind(mat, scores)
            }  ## for each gene section
        }  ## for each dataset

        browser()
        gg <- ggplot(melt(mat), aes(name="", x=X2, y=value))
        gg <- gg + geom_line(aes(color=X1))
        
        dimnames(mat) <- list(NULL, NULL)
        message(paste0("Dimensions of scoring matrix:", dim(mat)))
        d <- dist(mat) 
        message(paste0("Dimensions of square matrix:", dim(d)))
        fit <- cmdscale(d, eig=TRUE, k=2)
        x <- fit$points[,1]
        y <- fit$points[,2]

        gg <- ggplot(data.frame(x=x, y=y), aes(name="", x=x, y=y, label=datasets))
        gg <- gg + theme_bw()
        gg <- gg + theme(legend.position="top",legend.title=element_blank(), legend.text = element_text(colour="blue", size = 14, face = "bold"))
        gg <- gg + geom_point() 
        gg <- gg + geom_text()
        gg <- gg + scale_x_continuous(expand = c(.2,0))
        print(gg)
        
        ##plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",   main="Metric  MDS")
        ##text(x, y, labels = datasets, cex=.7)
    } ## for-each gene
}

main()
