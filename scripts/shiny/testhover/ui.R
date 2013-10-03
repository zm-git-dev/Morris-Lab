require(shiny)

shinyUI(
    pageWithSidebar(
        headerPanel("Test!"),
        sidebarPanel(
            checkboxInput(inputId = "noPlot",
                          label = "Change plot scale",
                          value = FALSE),
            p("Hover coordinates:"),
            strong(textOutput("coords")),
            helpText("Coordinates are always reported in the scale used by the first plot shown, regardless of how the plot is later modified.   Hover over a known point of the plot and note that the reported coordinate is correct.   Then check the box to change the scale of the graph.   Hover coordinates are *STILL* reported in the scale of the first plot!")
            ),
        mainPanel(
            plotOutput("plot",
                       hoverId="hover",
                       hoverDelay=800,
                       hoverDelayType="debounce")
            )
        )
    )

