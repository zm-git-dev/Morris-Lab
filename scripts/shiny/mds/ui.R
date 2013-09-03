

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
                            selected = "Human CLL"),
                p(strong("Select gene to analyze")),
                uiOutput('geneSelect')
                
                ## p(strong("Date range (months back from present);")),
                ## sliderInput(inputId = "obs",
                ##             label=" ",
                ##             min = 0, max = 60, step = 1, value = c(0,2))
                
                ), 
            
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
            tabsetPanel(
                tabPanel("Read Depth", plotOutput("rdplot")), 
                tabPanel("MDS", plotOutput("plot")), 
                tabPanel("Table", tableOutput("view"))
                )
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
