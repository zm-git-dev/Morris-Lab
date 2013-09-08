

library(shiny)


shinyUI(
    pageWithSidebar(
        
        ## Application title
        headerPanel("Characterizing ribosome profiles"),
        
        ## Sidebar with controls to select the subjects and time span
        sidebarPanel(

            helpText(p(
                "The graph represents number of aligned reads found on a particular gene in the specified datasets.")
                     ),
            
            wellPanel(
                p(strong("Select a set of ribomsome profiles to analyze")),
                selectInput(inputId = "dataspec",
                            label = "Profile data sets:",
                            choices = c(
                                "Mouse Prostate"="prostate",
                                "Human CLL"="cll",
                                "PEO1-RPT-corrUp"="PEO1-RPT-corrUp",
                                "PEO1-RPT-corrDown"="PEO1-RPT-corrDown",
                                "PEO1-RPT-top200"="PEO1-RPT-top200"
                                ),
                            selected = "PEO1-RPT-corrUp"),
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


            div(class="span6",
                actionButton("getGraph", "Get Graph")),
            div(class="span6", checkboxInput(inputId = "log", label = "log10 scale", value = FALSE)),
            
            helpText("Use log scale if compared searches are significantly different"),
            downloadButton('downloadData', 'Download Output as csv')
            
            ),

        
        ## Show the caption a line graph of the daily rate and summary of results 
        mainPanel(
            tags$head( tags$link(rel="stylesheet", type="text/css", href="app.css")),
            tabsetPanel(
                tabPanel("Read Depth",
                         div(id="plot_container",
                            plotOutput("rdplot"),
                             ## bootstrap dropdown menu code shamelessly stolen from
                             ## http://stackoverflow.com/a/13998987/1135316
                             ##
                             div(id="navicon", class="dropdown btn-group",
                                 a(img(src="navicon.svg"), class="btn btn-small dropdown-toggle", "data-toggle"="dropdown", href="#"),
                                 tags$ul(class="dropdown-menu pull-right",
                                         tags$li(a("Print chart", href="#", target="_blank")),
                                         tags$li(class="divider"),
                                         tags$li(a("Download PNG image", href="#", target="_blank")),
                                         tags$li(a("Download JPEG image", href="#", target="_blank")),
                                         tags$li(a("Download PDF Document", href="#", target="_blank"))
                                         )))),
                
                
                tabPanel("MDS", plotOutput("plot")), 
                tabPanel("Table", tableOutput("view"))
                )
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
