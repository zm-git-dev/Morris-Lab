require(shiny)

source("dropbutton.R")

shinyUI(
    pageWithSidebar(
        headerPanel("Test!"),
        sidebarPanel(
            list(div(HTML("I like <u>turtles</u>")),
                 textOutput("coords")),
            downloadButton('downloadData', 'Download Output as csv')
            ),
        mainPanel(
            h4(textOutput("caption")),
            tags$head( tags$link(rel="stylesheet", type="text/css", href="css/app.css")),
            tabsetPanel(
                tabPanel("ggplot",
                         with(tags, 
                              div(class="plot_container",
                                  plotOutput("radplot",
                                             hoverId="click",
                                             hoverDelay=800,
                                             hoverDelayType="debounce"),
                                  dropButton(inputId = "printmenu1",
                                             label = tags$img(src="navicon.svg"),
                                             choices = c(
                                                 "print"="Print chart",
                                                 "png"="Download PNG image",
                                                 "pdf"="Download PDF document"),
                                             class="pull-right navicon"
                                             
                                             )
                                  
                                  )
                              )
                         ),
                tabPanel("plot",
                         with(tags, 
                              div(class="plot_container",
                                  plotOutput("pplot",
                                             hoverId="click",
                                             hoverDelay=800,
                                             hoverDelayType="debounce"),
                                  dropButton(inputId = "printmenu2",
                                             label = tags$img(src="navicon.svg"),
                                             choices = c(
                                                 "print"="Print chart",
                                                 "png"="Download PNG image",
                                                 "pdf"="Download PDF document"),
                                             class="pull-right navicon"
                                             )
                                  )
                              )
                         )
                )
            )
        )
    )

