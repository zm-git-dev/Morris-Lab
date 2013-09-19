library(shiny)
library(ggplot2)
library(grid)

shinyServer(function(input, output, session) {
    output$caption <- renderText({
        paste0("printmenu: ", input$printmenu1)
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
        message(paste("plot range = [", paste(ggb$panel$ranges[[1]]$x.range, collapse=","), "]"))
        message(paste("plot range = [", paste(ggb$panel$ranges[[1]]$y.range, collapse=","), "]"))
    })

})
