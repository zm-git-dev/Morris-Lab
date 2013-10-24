

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
                                             "PEO1-RPT-top200"="PEO1-RPT-top200",
                                             "PEO1-RPT-corrUp"="PEO1-RPT-corrUp",
                                             "PEO1-RPT-corrDown"="PEO1-RPT-corrDown",
                                             )),
                             p(strong("Select gene to analyze")),
                             uiOutput('geneSelect')
                             
                             ## p(strong("Date range (months back from present);")),
                             ## sliderInput(inputId = "obs",
                             ##             label=" ",
                             ##             min = 0, max = 60, step = 1, value = c(0,2))
                             
                             ),

                         ## It would be nice to have a canvas covering the plot area with
                         ## a message and a busy indicator.
                         
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
                             checkboxInput(inputId = "colorOption",
                                           label = "Alternate palette",
                                           value = FALSE)),
                         helpText("Use a color-blind friendly palette."),

                         div(class="span6",
                             checkboxInput(inputId = "showSplices",
                                           label = "Show splice junctions",
                                           value = TRUE)),
                         helpText("Indicate splice junctions on read-depth plot."),

                         div(class="span6",
                             checkboxInput(inputId = "showLabels",
                                           label = "Label points",
                                           value = FALSE)),
                         helpText("Show labels on scatter plots"),

                         div(class="span6",
                             checkboxInput(inputId = "bindata",
                                           label = "Bin the data",
                                           value = TRUE)),
                         helpText("Speed up plotting by reducing plot fidelity"),

                         div(class="span6",
                             checkboxInput(inputId = "log",
                                           label = "Log scale",
                                           value = FALSE)),
                         helpText("Take log2 of read-depth data."),
                         
                         div(class="span6",
                             checkboxInput(inputId = "centerdata",
                                           label = "Mean center",
                                           value = FALSE)),
                         helpText("Center the read depth data around the mean"),

                         div(class="span6",
                             checkboxInput(inputId = "lowreads",
                                           label = "Filter low reads",
                                           value = FALSE)),
                         helpText("Only consider regions with minimum read depth"),

                         div(class="span6",
                             checkboxInput(inputId = "usenorm",
                                           label = "Pre-normalized data",
                                           value = TRUE)),
                         helpText("Use pre-normalized data if available"),

                         div(class="span6",
                             checkboxInput(inputId = "renorm",
                                           label = "Re-normalized data",
                                           value = FALSE)),
                         helpText("Divide each dataset by its sum")

                         )
                )
            ),
            
        
        ## Show tabbed panel with various graphs and tables.
        mainPanel(
            tags$head( tags$link(rel="stylesheet", type="text/css", href="css/app.css")),
            tabsetPanel(
                tabPanel("Read Depth",
                         plotOutput("rdplot"),
                         uiOutput('viewSlider'),
                         conditionalPanel(
                             condition="$('div#profile').hasClass('recalculating')",
                             img(src="loading.gif"))
                         ),
                tabPanel("MDS",
                         plotOutput("mdsplot", click="click")
                         ), 
                tabPanel("Table", tableOutput("view"))
                ) ## end-tabsetPanel
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
