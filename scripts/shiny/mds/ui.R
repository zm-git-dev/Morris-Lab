

require(shiny)
source("dropbutton.R")

shinyUI(
    pageWithSidebar(
        
        ## Application title
        headerPanel(list("Characterizing ribosome profiles",
                         img(src="FHlogo.gif", width=200, class="pull-right")),
                    windowTitle="Profiles"),
        
        ## Sidebar with controls to select the subjects and time span
        sidebarPanel(
            tabsetPanel(
                tabPanel("Datasets",
                         helpText(p(paste0(
                             "The graph represents number of aligned reads found on a particular",
                             "gene in the specified datasets."))),
                         wellPanel(
                             p(strong("Select a set of ribosome profiles to analyze")),
                             selectInput(inputId = "dataspec",
                                         label = "Profile data sets:",
                                         choices = c(
                                             "Mouse Prostate"="prostate",
                                             "Human CLL"="cll",
                                             "PEO1-RPT-corrUp"="PEO1-RPT-corrUp",
                                             "PEO1-RPT-corrDown"="PEO1-RPT-corrDown",
                                             "PEO1-RPT-top200"="PEO1-RPT-top200",
                                             "test"="test"
                                             ),
                                         selected = "test"),
                             p(strong("Select gene to analyze")),
                             uiOutput('geneSelect')
                             
                             ## p(strong("Date range (months back from present);")),
                             ## sliderInput(inputId = "obs",
                             ##             label=" ",
                             ##             min = 0, max = 60, step = 1, value = c(0,2))
                             
                             ),

                         ## It would be nice to have a canvas covering the plot area with
                         ## a message and a busy indicator.
                         
                         conditionalPanel(condition="$('div#rdplot').hasClass('recalculating')", img(src="loading.gif")),
                         conditionalPanel(condition="!($('div#rdplot').hasClass('recalculating'))", br()),


                         downloadButton('downloadData', 'Download Output as csv')
                         ),
                tabPanel("Options",
                         helpText("Select options to control the analysis of ribosome profiles."),

                         helpText("Distance algorithm for multi-dimensional sampling"),
                         selectInput(inputId = "distOption",
                                     label = "",
                                     choices = c(
                                         "euclidean", "maximum", "manhattan", "canberra",
                                         "binary", "minkowski")
                                     ),

                         p(),
                         div(class="span6",
                             checkboxInput(inputId = "log",
                                           label = "log scale",
                                           value = FALSE)),
                         helpText("Use log2 scale for read-depth plot."),
                         
                         div(class="span6",
                             checkboxInput(inputId = "colorOption",
                                           label = "Alternate palette",
                                           value = FALSE)),
                         helpText("Use a color-blind friendly palette."),

                         div(class="span6",
                             checkboxInput(inputId = "showSplices",
                                           label = "Show splice junctions",
                                           value = FALSE)),
                         helpText("Indicate splice junctions on read-depth plot."),

                         div(class="span6",
                             checkboxInput(inputId = "showLabels",
                                           label = "Label points",
                                           value = FALSE)),
                         helpText("Show labels on scatter plots"),

                         
                         helpText("Count read depth at a single point, or across the width of each read."),
                         selectInput(inputId = "rpsOption",
                                     label = "",
                                     choices = c(
                                         "Read depth"="stacked",
                                         "Centered point"="centered")
                                     )
                         )
                )
            ),
            
        
        ## Show tabbed panel with various graphs and tables.
        mainPanel(
            ## h3(textOutput("debug")),
            tags$head( tags$link(rel="stylesheet", type="text/css", href="css/app.css")),
            tabsetPanel(
                tabPanel("Read Depth",
                         with(tags, 
                              div(class="plot_container", plotOutput("rdplot"),
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
                
                
                tabPanel("MDS", 
                         with(tags, 
                              div(class="plot_container", plotOutput("mdsplot"),
                                  dropButton(inputId = "printmenu2",
                                             label = tags$img(src="navicon.svg"),
                                             choices = c(
                                                 "print"="Print chart",
                                                 "png"="Download PNG image",
                                                 "pdf"="Download PDF document")
                                             )
                                  )
                              )
                         ), 
                tabPanel("Table", tableOutput("view"))
                ) ## end-tabsetPanel
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
