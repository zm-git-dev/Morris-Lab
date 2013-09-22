# libraries used. install as necessary

# Time-stamp: <2013-09-21 00:43:27 chris>

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

        if (nrow(kg) != 1)
            return()

        ## remember the refseq name because that is what identifies each gene in a dataset
        refseq = kg[1,'name']

        print(paste0("retrieving data for ", dataset))
        df = morris.getalignments(dataset, refseq)
        attr(df, "dataset") <- descriptions[dataset,"description"]
        p = profile(df, kg[refseq,])

    })

    output$nrows <- reactive({
        d = data()
        return (if (is.null(d)) 0 else nrow(d))
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

        dataset = input$dataSelect
        if (is.null(dataset))
            return()
        stats = morris.fetchstats(dataset)

        df <- p$plotpositions()
        
        ## count how many reads occur on each position.
        histdata <- hist(df$rposition, breaks=c(1:p$transcript()$txLength()), plot=FALSE)
        par(mar=c(3, 3, 0.5, 1))  # Trim margin around plot [bottom, left, top, right]

        par(mgp=c(1.5, 0.2, 0))  # Set margin lines; default c(3, 1, 0) [title,labels,line]
        par(xaxs="r", yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

        if (input$normalize) {
            ## normalize the count at each position by the total RPM
            ## of mapped reads in the dataset.
            histdata$counts <- histdata$counts / (stats[dataset,"aligned_count"]/1e6)
        }
        xlim <- c(1,p$transcript()$txLength())
        
        plot(histdata$mids[histdata$counts != 0],histdata$counts[histdata$counts != 0],
             type='h', xlim=xlim, lwd=3, lend=2,
             xlab="", ylab="", cex=1.5, frame.plot=F, xaxt="n")
        print(par("yaxp"))
        print(paste0(c(0, max(histdata$counts), 6), collapse=", " ))
        ticks <- pretty(1:p$transcript()$txLength(), 4)
        print(paste0(ticks, collapse=", "))

        # par(mgp=c(axis.title.position, axis.label.position, axis.line.position))
        mgp <- par("mgp")
        print(paste0(mgp, collapse=", "))
        mgp[2] <- 0.5
        par(mgp=mgp)

        axis(1, at=ticks, labels=T, lwd=1, lwd.ticks=1, lty="solid", las=1, cex.axis=0.9)
        grid()

        if (input$showSplices) {
            ## Add an alpha value to a colour
            add.alpha <- function(col, alpha=1){
                if(missing(col))
                    stop("Please provide a vector of colours.")
                apply(sapply(col, col2rgb)/255, 2, 
                      function(x) 
                      rgb(x[1], x[2], x[3], alpha=alpha))  
            }

            lightred <- add.alpha("red", alpha=0.7)

            j <- p$transcript()$junctions()
            sapply(j, function(x) abline(v = x, col=lightred, lty="dashed"))
        }
    })
    
    ## draw the transcript anf the profile of reads along the transcript
    output$profile2 <- renderPlot({
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
