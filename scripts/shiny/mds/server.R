# libraries used. install as necessary

# Time-stamp: <2013-09-10 00:23:29 chris>

  library(shiny)
  library(ggplot2) # graphs
  library(plyr)  # manipulating data
  library(reshape)

## FIXME - this hardcoded path will not be portable.
  source("~/Morris-Lab/scripts/R/morrislib.R")
  source("~/Morris-Lab/scripts/R/profiling/profiling.R")


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

shinyServer(function(input, output, session) {
  


    ## Reactive function to create a list of datasets.  This list may
    ## contain a set of specifications for a database, or it might be a
    ## specification for retrieving data from a file.
    inputDatasets <- reactive({
        if (is.null(input$dataspec))
            return()
        print("entering data reactive function")
        ## Choose the datasets that meet the user's criteria
        if (input$dataspec == "prostate") {
            control <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+")
            treated <- c("103112_A", "030713_B", "041713_B")
            datasets = c(control=control, treated=treated)
        } else if (input$dataspec == "cll") {
            control = c("061113_A")
            treated = c("061113_B", "061113_C", "061113_D")
            datasets = c(control=control, treated=treated)
        } else if (input$dataspec == "PEO1-RPT-corrUp") {
            control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
            treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "~/Downloads/PEO1-RPT-corrUp-coverage.tsv"
        } else if (input$dataspec == "PEO1-RPT-corrDown") {
            control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
            treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "~/Downloads/PEO1-RPT-corrDown-coverage.tsv"
        } else if (input$dataspec == "PEO1-RPT-top200") {
            control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
            treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "~/Downloads/PEO1-RPT-top200-coverage.tsv"
        } else if (input$dataspec == "test") {
            control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
            treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "~/Downloads/test.tsv"
        } else {
            stop(paste0("unknown data specification - ", input$dataspec))
        }
        return(datasets)
    })


    ## Reactive function to retrieving complete data for a dataset
    dataFromFile <- reactive({
        ## only works for datasets with filename attribute
        datasets = inputDatasets()
        if (is.null(attr(datasets, "filename")))
            return()
        data = read.csv(attr(datasets, "filename"), sep='\t')
        data = data[, c("symbol", datasets)]
    })
        

    # Reactive function for retrieving a list of genes covered in the
    # datasets as specified by inputDatasets.
    geneChoices <- reactive({
        print("in gene choices")
        datasets = inputDatasets()
        if (is.null(datasets))
            return()

        # if the data is stored in a file, we'll pretty much have to read the whole thing in.
        if (! is.null(attr(datasets, "filename"))) {
            data = dataFromFile()
            ## extract the list of genes.
            return(levels(data$symbol))
        }

        ## Otherwise the data is stored in a database and we can let
        ## the DB filter out the gene list for us.
        
        ## Get the top expressed genes for the aggregate datasets?
        ## Get all the known genes for the genome?

        ## Try getting the top expressed genes.
        ## Take the datasets, grab the gene list from the dataset
        ## run unique on it.
        ## as.data.frame(table(datasets))
        else if (input$dataspec == "prostate") {
            return(c(c( "Tubb5", "Hist2h2bb", "Nov", "Chgb", "Hist1h1a" ),
                     c("Tgm4", "Pbsn", "Sbp", "Rnase1"),
                     c("Vim", "Itgb1", "Itga1", "Col1a1", "Col1a2")))
        } else  if (input$dataspec == "cll") {
            return(c("RPS2", "RPS23", "EIF3E", "EIF3H", "EIF3I", "EEF1a1", "EIF4EBP2"))
        } else {
            stop(paste0("unknown data specificartion - ", input$dataspec))
        }
        
    })
                
                
    output$geneSelect <- renderUI({
        print("in renderUI of geneselect")
        choices = geneChoices()
        if(is.null(choices))
            return()
        selectInput("geneSelect", "", choices=choices)
    })

    ## Reactive function for retrieving a matrix representing aligned
    ## read positions on a patricular gene for a collection of
    ## datasets.  The datasets might repesent replicated control and
    ## treated conditions, and the gene is one particular gene chosen
    ## from among those in the datasets
    ##
    ## FIXME: This reactive should be prioritized below that of
    ## geneChoices, as this routine cannot proceed until a valid gene
    ## has been selected from the list, and that list is dependent on
    ## the dataset specification.
    data <- reactive({
        gene = input$geneSelect
        if(is.null(gene))
            return()

        datasets = inputDatasets()
        if (is.null(datasets))
            return()
        
        if (! is.null(attr(datasets, "filename"))) {
            data <- dataFromFile()
            data <- data[data$symbol==gene,]
            rownames(data) = 1:nrow(data)
            mat = as.matrix(t(subset(data,,-symbol)))
            ## extract the list of genes.
            return(mat)
        }
        
        genome <- morris.getGenome(datasets[[1]])
        descriptions = morris.fetchinfo(datasets)[,"description", drop=FALSE]
        stats = morris.fetchstats(datasets)

        kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
        rownames(kg) <- kg$name
        stopifnot(nrow(kg) == 1)

        ## remember the refseq name because that is what identifies each gene in a dataset
        refseq = kg[1,'name']

        ## create an empty matrix to hold the experimental data
        mat <- matrix(0, 0, 50)
        for (dataset in datasets) {
            print(paste0("retrieving data for ", dataset))
            df = morris.getalignments(dataset, refseq)
            attr(df, "dataset") <- descriptions[dataset,"description"]
            ribo.profile = profile(df, kg[refseq,])
            df <- ribo.profile$plotpositions()
            xlim <- c(1, ribo.profile$transcript()$txLength())

            ## count how many reads occur on each position.
            histdata <- hist(df$rposition, breaks=c(1:ribo.profile$transcript()$txLength()), plot=FALSE)
            ## normalize the count at each position by the total RPM
            ## of mapped reads in the dataset.
            scores <- histdata$counts / (stats[dataset,"aligned_count"]/1e6)
            print(paste0("normalization factor for ", dataset, " = ", (stats[dataset,"aligned_count"]/1e6)))
            dim(mat) = c(dim(mat)[1], length(scores))
            mat <- rbind(mat, scores)
            print(paste0("finished data for ", dataset))
        }  ## for each dataset
        rownames(mat) = datasets
        
        print("returning data")
        return(mat)
    })
    
    
    ## Create a heading based on range of dates selected for printing as a caption
    output$caption <- renderText({
        foo = input$printmenu1
        paste("", Sys.time())
    })
    
    
    ## create 2-D scatter plot of multidimensional sampling data
    output$mdsplot <- renderPlot({
        mat = data()
        if (is.null(mat))
            return()

        ## get the datasets, as this will tell us which are control
        ## group and which are the experimental group.
        datasets = inputDatasets()
        if (is.null(datasets))
            return()
        treated <- grep("treated", names(datasets))
        control <- setdiff(1:length(datasets), treated)
        
        condition <- sapply(names(datasets), function(x) !is.na(pmatch('treated',x)))

        d <- dist(mat, method=input$distOption)
        fit <- cmdscale(d, eig=TRUE, k=2)

        ## user can select a standard palette or one that is more
        ## visible to those with R-G color blinkdness.
        colorPalette <- switch(input$colorOption,
                               "std"=c("red", "blue"),
                               "rg-cb"=cbPalette)

        gg.data = data.frame("x"=fit$points[,1], "y"=fit$points[,2], color=condition)

        gg <- ggplot(gg.data, aes(name="", x=x, y=y, label=rownames(mat)),
                     environment = environment())
        gg <- gg + theme_bw()
        gg <- gg + theme(axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.title=element_blank(),
                         ## legend.position="bottom",
                         legend.direction="horizontal",
                         legend.text = element_text(colour="blue", size = 14, face = "bold"))
        gg <- gg + theme(legend.justification=c(.5,1), legend.position=c(.5,0))
        gg <- gg + theme(legend.key = element_rect(size = 0.5, linetype="blank"))
        gg <- gg + theme(panel.border = element_blank())
        ##  gg <- gg + theme(legend.background = element_rect(colour = 'purple', fill = 'pink', size = 1))
        
        gg <- gg + geom_point(size=4, aes(group=condition, color=condition)) 
        ##gg <- gg + geom_text()
        gg <- gg + scale_colour_manual(values=colorPalette,
                                       labels=c("Control", "Treated"))
        gg <- gg + scale_x_continuous(expand = c(.2,0))
        print(gg)
    })
    
    

    ## create line plot for read depth
    output$rdplot <- renderPlot({
        mat = data()
        if (is.null(mat))
            return()

        ## get the datasets, as this will tell us which are control
        ## group and which are the experimental group.
        datasets = inputDatasets()
        if (is.null(datasets))
            return()
        treated = grep("treated", names(datasets))
        control = setdiff(1:length(datasets), treated)
        
        gg <- ggplot(melt(mat), aes(name="", x=X2, y=value, group=X1))
        gg <- gg + theme_bw()
        ## gg <- gg + theme(legend.position="top",
        ##                  legend.title = element_text(colour="black", size = 14, face = "bold"),
        ##                  legend.text = element_text(colour="blue", size = 12, face = "bold"))
        ## gg <- gg + scale_colour_discrete(name = paste(input$dataspec, ":", input$geneSelect))
        gg <- gg + ylab("read depth (normalized to RPM)")
        gg <- gg + xlab("Transcript position")
        ##gg <- gg + geom_line(aes(color=X1))
        gg <- gg + geom_line(data=melt(mat[treated,,drop=FALSE]), aes(color="red"))
        gg <- gg + geom_line(data=melt(-mat[control,,drop=FALSE]), aes(color="blue"))
                                             "PEO1-RPT-top200"="PEO1-RPT-top200"
                                             "PEO1-RPT-top200"="PEO1-RPT-top200"
        print(gg)
    })
    
    
    
    ## create summary data for each subject 
    output$view <- renderTable({
        mat = data()
        if (is.null(mat))
            return()
        mat
    })
    
    ## make data downloadable
    output$downloadData <- downloadHandler(
        ## filename = function() { paste(input$data, '.csv', sep='') }, in tutorial to distinguish files trickier with my work
        filename = function() { paste('results.csv', sep='') },
        content = function(file) {
            write.csv(data(), file)
        }
        )
    
})
