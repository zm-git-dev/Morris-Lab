# libraries used. install as necessary

# Time-stamp: <2013-09-24 17:24:02 chris>

  library(shiny)
  library(plyr)  # manipulating data
  library(reshape)

## FIXME - this hardcoded path will not be portable.
  source("../shared/R/morrislib.R")
  source("../shared/R/profiling.R")


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

shinyServer(function(input, output, session) {

    observe({
        message("Inside printmenu observer")
        selection = input$printmenu1
        if (is.null(selection))
            return()
        selection = sub("^([^.]*)-.*", "\\1", selection)
        if (selection == "print") {
            message("	invoke print")
        } else if (selection == "png") {
            message("	invoke png")
        } else if (selection == "pdf") {
            message("	invoke pdf")
        } else if (selection == "stop") {
            message("	invoke stop")
            shiny::stopApp()
        } else {
            message(paste0("	unknown printmenu command - ", selection) )
        }
    })

    output$downloadData <- downloadHandler(
        filename = function() {
            paste('data-', Sys.Date(), '.png', sep='')
        },
        content = function(file) {
            png(file, width = 980, height = 400, units = "px", pointsize = 12,
                bg = "white", res = NA)
            output$profile()
            dev.off()
        },
        contentType = 'image/png'
        )



    # Reactive function for retrieving a list of genes covered in the
    # datasets as specified by inputDatasets.
    geneChoices <- reactive({
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

    rptData <- reactive({
        gene <- input$geneName
        if(is.null(gene))
            return()

        dataset <- input$dataSelect
        if (is.null(dataset)) 
            return()

        genome <- morris.getGenome(dataset)
        descriptions = morris.fetchinfo(dataset)[,"description", drop=FALSE]
        kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
        rownames(kg) <- kg$name

        if (nrow(kg) != 1)
            return()

        ## remember the refseq name because that is what identifies each gene in a dataset
        refseq <- kg[1,'name']

        print(paste0("retrieving data for ", dataset))
        df <- morris.getalignments(dataset, refseq)
        attr(df, "dataset") <- descriptions[dataset,"description"]
        p <- profile(df, kg[refseq,])

    })

    ## Create a heading for debugging....
    output$caption <- renderText({
        hover = input$plotHover
        sprintf("x=%.4g  y=%.4g", hover$x, hover$y)
    })

    ## draw the profile of reads along the transcript
    output$profile <- renderPlot({
        p = rptData()
        dataset = input$dataSelect
        if (is.null(p) || is.null(dataset)) 
            return()

        view = input$viewSlider
        if (is.null(view)) {
            warning("viewSlider returned NULL!")
            view <- c(1,p$transcript()$txLength())
        }

        df <- p$plotpositions()
        
        ## count how many reads occur on each position.
        histdata <- hist(df$rposition, breaks=c(1:p$transcript()$txLength()), plot=FALSE)
        par(mar=c(3, 3, 0.5, 1))  # Trim margin around plot [bottom, left, top, right]

        par(mgp=c(1.5, 0.2, 0))  # Set margin lines; default c(3, 1, 0) [title,labels,line]
        par(xaxs="i", yaxs="i")  # Extend axis limits by 4% ("i" does no extension)

        if (input$normalize) {
            ## normalize the count at each position by the total RPM
            ## of mapped reads in the dataset.
            stats = morris.fetchstats(dataset)
            histdata$counts <- histdata$counts / (stats[dataset,"aligned_count"]/1e6)
        }
        xlim <- view

        
        plot(histdata$mids[histdata$counts != 0],histdata$counts[histdata$counts != 0],
             type='h', xlim=xlim, lwd=3, lend=2,
             xlab="", ylab="", cex=1.5, frame.plot=F,
             xaxt="n"
             )

        # Draw a pretty X-axis.  The default axis sometimes does
        # screwy things, like not spanning the entire width of the
        # plot.  By drawing our own we have more control over exacty
        # how it looks.
        #
        # par(mgp=c(axis.title.position, axis.label.position, axis.line.position))
        ticks <- pretty(view, 4)
        mgp <- par("mgp")
        mgp[2] <- 0.5
        par(mgp=mgp)
        axis(1, at=ticks, labels=T, lwd=1, lwd.ticks=1, lty="solid", las=1, cex.axis=0.9)

        # draw a subtle grid in the background of the plot
        grid(nx=0)

        # show vertical lines at splice sites of the transcript.
        # user can turn this off on the 'options' tab.
        if (input$showSplices) {
            # Add an alpha value to a colour
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

        # annotate the graph with the name of the dataset and the number of reads.
        # The usr coordinates appear to change after plot(), so fetch the coordinates
        # again before positioning the labels.
        tmp = par("usr")
        text(xlim[1], tmp[4], attr(df, "dataset"), adj=c(0,1.2), new=TRUE)
        text(xlim[2], tmp[4], paste0(nrow(df), " reads"), adj = c( 1, 1.2 ), new=TRUE)

    })
    
    output$viewSlider <- renderUI({
        p = rptData()
        if (is.null(p) || class(p) != "profile")
            return()

        sliderInput(inputId = "viewSlider",
                    label=" ",
                    min = 1, max = p$transcript()$txLength(), step = 1,
                    value = c(1,p$transcript()$txLength()))
    })

    ## create summary data for each subject 
    output$view <- renderTable({
        mat = rptData()
        if (is.null(mat))
            return()
        mat
    })
    
    ## make data downloadable
    output$downloadData <- downloadHandler(
        ## filename = function() { paste(input$data, '.csv', sep='') }, in tutorial to distinguish files trickier with my work
        filename = function() { paste('results.csv', sep='') },
        content = function(file) {
            write.csv(rptData(), file)
        }
        )
    
})
