# libraries used. install as necessary

# Time-stamp: <2013-08-28 13:39:41 chris>

  library(shiny)
  library(RJSONIO) # acquiring and parsing data
  library(ggplot2) # graphs
  library(plyr)  # manipulating data
  library(lubridate) #dates
  library(stringr)
  library(reshape)

## FIXME - this hardcoded path will not be portable.
  source("~/Morris-Lab/scripts/R/morrislib.R")
  source("~/Morris-Lab/scripts/R/profiling/profiling.R")

trim.leading <- function (x)  sub("^\\s+", "", x)

shinyServer(function(input, output, session) {
  
    geneChoices <- reactive({
        print("in gene choices")
        if(is.null(input$dataspec))
            return()

        ## Get the top expressed genes for the aggregate datasets?
        ## Get all the known genes for the genome?

        ## Try getting the top expressed genes.
        ## Take the datasets, grab the gene list from the dataset
        ## run unique on it.
        ## as.data.frame(table(datasets))
        if (input$dataspec == "prostate") {
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

    inputDatasets <- reactive({
        print("entering data reactive function")
        ## Choose the datasets that meet the user's criteria
        if (input$dataspec == "prostate") {
            control <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+")
            treated <- c("103112_A", "030713_B", "041713_B")
        } else if (input$dataspec == "cll") {
            control = c("061113_A")
            treated = c("061113_B", "061113_C", "061113_D")
        } else {
            stop(paste0("unknown data specificartion - ", input$dataspec))
        }
        c(treated,control)
    })

    ## a reactive conductor function, carries out a long-running computation.
    data <- reactive({
        if(is.null(input$geneSelect))
            return()

        datasets = inputDatasets()
        genome <- morris.getGenome(datasets[[1]])
        descriptions = morris.fetchinfo(datasets)[,"description", drop=FALSE]
        stats = morris.fetchstats(datasets)

        gene = input$geneSelect
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
        if(is.null(input$geneSelect) || is.null(input$dataspec))
            return()
        paste(input$dataspec, ":", input$geneSelect)
    })
    
    
    ## create 2-D scatter plot of multidimensional sampling data
    output$plot <- renderPlot({
        mat = data()
        if (is.null(mat))
            return()
        d <- dist(mat)
        print("calling cmdscale")
        fit <- cmdscale(d, eig=TRUE, k=2)
        x <- fit$points[,1]
        y <- fit$points[,2]
        gg <- ggplot(data.frame(x=x, y=y), aes(name="", x=x, y=y, label=rownames(mat)),
                     environment = environment())
        gg <- gg + theme_bw()
        gg <- gg + theme(legend.position="top",legend.title=element_blank(), legend.text = element_text(colour="blue", size = 14, face = "bold"))
        gg <- gg + geom_point() 
        gg <- gg + geom_text()
        gg <- gg + scale_x_continuous(expand = c(.2,0))
        print(gg)
    })
    
    

    ## create line plot for read depth
    output$rdplot <- renderPlot({
        mat = data()
        if (is.null(mat))
            return()
        gg <- ggplot(melt(mat), aes(name="", x=X2, y=value))
        gg <- gg + theme_bw()
        gg <- gg + theme(legend.position="right",
                         legend.title = element_text(colour="black", size = 14, face = "bold"),
                         legend.text = element_text(colour="blue", size = 12, face = "bold"))
        gg <- gg + scale_colour_discrete(name = paste(input$dataspec, ":", input$geneSelect))
        gg <- gg + ylab("read depth (normalized to RPM)")
        gg <- gg + xlab("Transcript position")
        gg <- gg + geom_line(aes(color=X1))
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
