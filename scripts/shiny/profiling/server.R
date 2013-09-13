# libraries used. install as necessary

# Time-stamp: <2013-09-12 14:35:55 chris>

  library(shiny)
  library(plyr)  # manipulating data
  library(reshape)

## FIXME - this hardcoded path will not be portable.
  source("../shared/R/morrislib.R")
  source("../shared/R/profiling.R")


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

shinyServer(function(input, output, session) {


    # Reactive function for retrieving a list of genes covered in the
    # datasets as specified by inputDatasets.
    geneChoices <- reactive({
        print("in gene choices")
        dataset = input$dataSelect
        if (is.null(dataset))
            return()

        genome <- morris.getGenome(dataset)
        df = morris.genecounts(dataset)

        topGenes = rownames(df[order(df[,1], decreasing=TRUE)[1:10],,drop=FALSE])
        kg = morris.getknowngenes(genome, gene=topGenes)
        genes = kg$name2
        
    })
                
                
    output$orgSelect <- renderUI({
        datasets = morris.datasets()
        if(is.null(datasets))
            return()
        df = morris.fetchinfo(datasets)
        organisms = unique(df$organism)
        selectInput(inputId = "orgSelect",
                    label="",
                    choices=organisms)
    })

    output$dataSelect <- renderUI({
        org = input$orgSelect
        if(is.null(org))
            return()
        datasets = morris.datasets(organism=org)
        if(is.null(datasets))
            return()
        descriptions = morris.fetchinfo(datasets)[,"description", drop=FALSE]
        descriptions = descriptions[!is.null(descriptions)]

        ## convert single-column data frame to a vector with names
        choices = rownames(descriptions)
        names(choices) = descriptions$description
        ## create the data selection control
        selectInput(inputId = "dataSelect",
                    label = "Profile data sets:",
                    choices=choices)
    })


    output$geneSelect <- renderUI({
        print("in renderUI of geneselect")
        choices = geneChoices()
        if(is.null(choices))
            return()
        selectInput("geneSelect", "", choices=choices)
    })

    data <- reactive({
        gene = input$geneSelect
        if(is.null(gene))
            return()

        dataset = input$dataSelect
        if (is.null(dataset))
            return()

        genome <- morris.getGenome(dataset)
        descriptions = morris.fetchinfo(dataset)[,"description", drop=FALSE]
        kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
        rownames(kg) <- kg$name
        stopifnot(nrow(kg) == 1)

        ## remember the refseq name because that is what identifies each gene in a dataset
        refseq = kg[1,'name']

        print(paste0("retrieving data for ", dataset))
        df = morris.getalignments(dataset, refseq)
        attr(df, "dataset") <- descriptions[dataset,"description"]
        p = profile(df, kg[refseq,])

    })
    
    
    ## Create a heading based on range of dates selected for printing as a caption
    output$caption <- renderText({
        if(is.null(input$geneSelect) || is.null(input$dataSelect))
            return()
        paste(input$dataSelect, ":", input$geneSelect)
    })
    
    ## draw the transcript anf the profile of reads along the transcript
    output$profile <- renderPlot({
        p = data()
        if (is.null(p))
            return()

        rpos = plot(p, minlen=0, units="aa")
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
