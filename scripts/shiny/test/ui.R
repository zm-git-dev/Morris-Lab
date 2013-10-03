require(shiny)

source("dropbutton.R")

shinyUI(
    pageWithSidebar(
        headerPanel("Test!"),
        sidebarPanel(
            list(div(HTML("I like <u>turtles</u>")))
            ),
        mainPanel(
            h4(textOutput("caption")),
            tags$head( tags$link(rel="stylesheet", type="text/css", href="css/app.css")),
            tabsetPanel(
                tabPanel("Radiation Exposure",
                         with(tags, 
                              div(class="plot_container", plotOutput("radplot"),
                                  dropButton(inputId = "printmenu1",
                                             label = tags$img(src="navicon.svg"),
                                             choices = c(
                                                 "print"="Print chart",
                                                 "png"="Download PNG image",
                                                 "pdf"="Download PDF document")
                                             )
                                  )
                              )
                         )
                )
            )
        )
    )

