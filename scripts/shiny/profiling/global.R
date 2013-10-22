

getCDSAlignments <- function(dataset, group=NULL) {
    ## this routine operates only on ONE dataset name - NOT a vector!
    stopifnot(length(dataset) == 1)

    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop("Could not connect to database: ", e$message);
        }

        ## check to see that each experiment has the same genome!
        ## it won't be useful to compare datasets aligned with different genomes. 
        ## ( actually this routine only works with a single dataset, so really there is no
        ##   need to consolidate the results from multiple datasets. )
        gs <- morris.getGenome(dataset)
        genome <- unique(gs)
        if (length(genome) == 0) {
            stop("Error during database query: no dataset found with the name ", paste0(dataset, collapse=", "))
        } else if (length(genome) != 1) {
            stop("Error during database query: not all selected datasets align to the same genome")
        }

        ## SQL query for retrieving the per-gene alignment count of a single dataset.
        ## These individual tables will be joined together based on the gene name
        ## as a common key.   Oddly eough it is faster to do this in R than it is
        ## in SQL.

        query <- paste(sprintf('select a.feature, k.name2 as "common", count(*) as "%s" from ', dataset),
                       '( morris.new_alignments_tbl a join datasets_tbl d ON a.dataset_id = d.id)',
                       'join knowngenes2 k on k.name=a.feature',
                       sprintf('where d.name like "%s%%"', dataset),
                       'and (a.position+(case when (a.strand = "+") then (a.length) else (-a.length) end)/2)',
                       'between k.cdsStart and k.cdsEnd',
                       'group by a.feature')
        df <- dbGetQuery(con, query)
        if (nrow(df) == 0) {
            stop("Error during database query: no alignments in dataset '", dataset,
                 "' for gene(s) ", paste0(gene, collapse=","), "\n",
                 "Available datasets:\n",
                 as.character(morris.datasets(group)), call.=TRUE)
        }
        rownames(df) <- df$feature
        df <- subset(df, select=c(-feature))
        attr(df, "genome") <- genome
        attr(df, "dataset") <- dataset
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })

    return(df)


}
