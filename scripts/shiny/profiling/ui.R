

library(shiny)
source("../shared/R/dropbutton.R")


shinyUI(
    pageWithSidebar(
        
        ## Application title
        headerPanel(list("Profiles", img(src="UWlogo.jpg", width=200, class="pull-right")),
                    windowTitle="Profiles"),
        
        ## Sidebar with controls to select the subjects and time span
        sidebarPanel(
            tags$head( tags$link(rel="stylesheet", type="text/css", href="css/app.css")),
            tags$head( tags$script(src="app.js")),
            tabsetPanel(
                tabPanel("Datasets",
                         helpText(p(paste0(
                             "The graph represents number of aligned reads found on a particular",
                             "gene in the specified datasets."))),
                         wellPanel(
                             uiOutput('orgSelect'),
                             uiOutput('dataSelect'),
                             p(strong("Select gene to analyze")),
                             uiOutput('geneSelect')
                             ),

                         ## It would be nice to have a canvas covering the plot area with
                         ## a message and a busy indicator.
                         
                         conditionalPanel(condition="$('div#rdplot').hasClass('recalculating')", img(src="loading.gif")),
                         conditionalPanel(condition="!($('div#rdplot').hasClass('recalculating'))", br()),

                         downloadButton('downloadData', 'Download Output as csv')

                         ),
                tabPanel("Options",
                         helpText((paste0(
                             "Select options to control the analysis of ribosome profiles."))),
                         wellPanel(
                             div(class="span6",
                                 checkboxInput(inputId = "log",
                                               label = "Log scale",
                                               value = FALSE)),
                             helpText("Use log2 scale for read-depth plot."),
                             
                             div(class="span6",
                                 checkboxInput(inputId = "colorOption",
                                               label = "Alternate palette",
                                               value = FALSE)),
                             helpText("Use a color-blind friendly palette."),

                             div(class="span6",
                                 checkboxInput(inputId = "exonBoundaries",
                                               label = "Show boundaries",
                                               value = FALSE)),
                             helpText("Show exon boundaries on transcript graph.")


                             )
                         )

                )
            
            ),
        
        ## Show the caption a line graph of the daily rate and summary of results 
        mainPanel(
            tabsetPanel(
                tabPanel("Read Depth",
                         with(tags, 
                              div(class="plot_container", plotOutput("profile"),
                                  dropButton(inputId = "printmenu1",
                                             label = tags$img(src="navicon.svg"),
                                             choices = c(
                                                 "print"="Print chart",
                                                 "png"="Download PNG image",
                                                 "pdf"="Download PDF document")
                                             )
                                  )
                              )
                         ),
                tabPanel("Table", tableOutput("view")),
                tabPanel("NEWS",
                         includeHTML("NEWS.html"))
                )
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
