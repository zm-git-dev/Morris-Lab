library(shiny)

shinyServer(function(input, output, session) {
    output$coords = renderText({
        hover = input$hover
        sprintf("x=%.4g   y=%.4g", hover$x, hover$y)
    })

    output$plot <- renderPlot({ 
        cars <- c(1, 3, 6, 4, 9)
        if (input$noPlot) {
            cars <- cars * 100
        }
        plot(cars)
    })
})
