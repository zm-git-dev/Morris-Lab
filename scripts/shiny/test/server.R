library(shiny)
library(ggplot2)
library(grid)

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

    output$caption <- renderText({
        paste0("printmenu: ", input$printmenu1)
    })

    output$coords = renderText({
        sprintf("hi %.4g %.4g", (input$click)$x, (input$click)$y)

    })

    output$radplot <- renderPlot({
        y <- rt(200, df = 5)
        gg = qplot(sample = y, stat="qq")
        gg <- gg + theme_bw()
        gg <- gg + theme(legend.title = element_text(size = 16, face = "bold"),
                         legend.text = element_text(size = 14, face = "bold"),
                         legend.position="top",
                             legend.direction="horizontal")
        gg <- gg + theme(legend.key = element_rect(size = 0.5, linetype="blank"))
        gg <- gg + theme(panel.border = element_blank())
        gg <- gg + theme(plot.margin = unit(c(0,0,0,0), "cm"))


        print(gg)
        ggb = ggplot_build(gg)
        message(paste("plot range X = [", paste(ggb$panel$ranges[[1]]$x.range, collapse=","), "]"))
        message(paste("plot range Y = [", paste(ggb$panel$ranges[[1]]$y.range, collapse=","), "]"))

    })

    output$pplot <- renderPlot({
        # Define the cars vector with 5 values
        cars <- c(1, 3, 6, 4, 9)

        view=c(-4,12)
        
        par(mar=c(3, 3, 0.5, 1))  # Trim margin around plot [bottom, left, top, right]

        par(mgp=c(1.5, 0.2, 0))  # Set margin lines; default c(3, 1, 0) [title,labels,line]
        par(xaxs="r", yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

        # Graph the cars vector with all defaults
        plot(cars, type='h',  lwd=3, lend=2,
             xlab="", ylab="", cex=1.5, frame.plot=F,
             xaxt="n")

        ticks <- pretty(view, 4)
        mgp <- par("mgp")
        print(paste0(mgp, collapse=", "))
        mgp[2] <- 0.5
        par(mgp=mgp)
        axis(1, at=ticks, labels=T, lwd=1, lwd.ticks=1, lty="solid", las=1, cex.axis=0.9)

        # draw a subtle grid in the background of the plot
        grid()



    })

        ## make data downloadable
    output$"menu-printmenu1-1" <- downloadHandler(
        ## filename = function() { paste(input$data, '.csv', sep='') }, in tutorial to distinguish files trickier with my work
        filename = function() { paste('results.csv', sep='') },
        content = function(file) {
            write.csv(data(), file)
        }
        )

        ## make data downloadable
    output$downloadData <- downloadHandler(
        filename = function() { message("in downloadHandler for downloadData button")
                                browser()
                                paste('results.csv', sep='') },
        content = function(file) {
            write.csv(data(), file)
        }
        )

})
