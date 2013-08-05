source("../morrislib.R")
source("profiling.R")

datasets= c("061113_A", "061113_B", "061113_C", "061113_D")
genome <- morris.getGenome(dataset[1])
descriptions = morris.fetchdesc(datasets)
motiflen = 9

profiles <- list(
    RPS23 = list(c())
)

score <- function(s, pssm) {
    ## For each character at each position in the sequence, add up the
    ## log-likelihood ratios found in the corresponding row and column
    ## of the pssm.
    return (sum(mapply(function(x,y) pssm[x,y], unlist(strsplit(s,'')), seq(nchar(s)))))
}

main <- function() {

    for (gene in names(profiles)) {
        kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
        rownames(kg) <- kg$name
        stopifnot(nrow(kg) == 1)
        
        refseq = kg[1,'name']

        for (dataset in datasets) {
            df = morris.getalignments(dataset, refseq)
            attr(df, "dataset") <- descriptions[dataset,"description"]
            ribo.profile = profile(df, kg[refseq,])
            print(paste0(refseq,":", ribo.profile$transcript()$name2()))
            df <- ribo.profile$plotpositions()
            print(head(df))
            browser()
            for (x in profiles[[gene]]) {
                if (is.null(x)) {
                    xlim <- c(1, ribo.profile$transcript()$txLength())
                } else {
                    xlim = as.numeric(x)
                }
                
                fscore =  function(i) {
                    with(df, sum( (rpositions >= i) 
                                  & (rpositions <= (i+motiflen-1))))
                     }
                scores = sapply(1:( ribo.profile$transcript()$txLength()), fscore)
                cat(scores)
            }  ## for each gene section
        }  ## for each dataset
    } ## for-each gene
}

main()
