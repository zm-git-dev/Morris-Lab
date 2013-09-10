

library(shiny)



shinyUI(
    pageWithSidebar(
        
        ## Application title
        headerPanel(list(img(src="UWlogo.jpg", width=200), "Profiles"),
                    windowTitle="Profiles"),
        
        ## Sidebar with controls to select the subjects and time span
        sidebarPanel(

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


                         div(class="span6", checkboxInput(inputId = "log",
                                 label = "log10 scale", value = FALSE)),
                         
                         helpText("Use log scale if read depth has significant vairance"),
                         downloadButton('downloadData', 'Download Output as csv')

                         ),
                tabPanel("Options",
                         helpText((paste0(
                             "Select options to control the analysis of ribosome profiles."))),
                         wellPanel(
                             selectInput(inputId = "colorOption",
                                         label = "Color pallete:",
                                         choices = c(
                                             "standard"="std",
                                             "R-G colorblind"="rg-cb")
                                         ),
                             selectInput(inputId = "rpsOption",
                                         label = "Read position:",
                                         choices = c(
                                             "Stacked reads"="stacked",
                                             "Centered reads"="centered")
                                         )
                             )
                         )

                )
            
            ),
        
        ## Show the caption a line graph of the daily rate and summary of results 
        mainPanel(
            tabsetPanel(
                tabPanel("Read Depth", plotOutput("profile")),
                tabPanel("Table", tableOutput("view"))
                )
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
