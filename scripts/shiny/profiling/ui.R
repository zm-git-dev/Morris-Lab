

library(shiny)
source("../shared/R/dropbutton.R")


shinyUI(
    pageWithSidebar(
        
        ## Application title
        headerPanel(list("Profiles", img(src="UWlogo.jpg", width=200, class="pull-right")),
                    windowTitle="Profiles"),
        
        ## Sidebar with controls to select the subjects and time span
        sidebarPanel(
            tags$head(
                tags$link(rel="stylesheet", type="text/css", href="css/app.css")),
            tabsetPanel(
                tabPanel("Datasets",
                         helpText(p(paste0(
                             "The graph represents number of aligned reads found on a particular",
                             "gene in the specified datasets."))),
                         wellPanel(
                             uiOutput("orgSelect"),
                             uiOutput('dataSelect'),
                             helpText(strong("Select gene to analyze"),
                                      strong("Enter gene name, correctly spelt")),
                             textInput(inputId = "geneName", label = " ", value = "Pbsn"),
                             downloadButton('downloadGene', "Download this gene")

                             ),

                         ## It would be nice to have a canvas covering the plot area with
                         ## a message and a busy indicator.
                         
                         conditionalPanel(condition="$('div#rdplot').hasClass('recalculating')", img(src="loading.gif")),
                         conditionalPanel(condition="!($('div#rdplot').hasClass('recalculating'))", br()),
                         downloadButton('downloadDataset', "Download all")

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
                                 checkboxInput(inputId = "showSplices",
                                               label = "Show junctions",
                                               value = TRUE)),
                             helpText("Show exon boundaries on transcript graph."),

                             div(class="span6",
                                 checkboxInput(inputId = "useCodons",
                                               label = "Use codons",
                                               value = FALSE)),
                             helpText("Aggregate reads into codons."),

                             div(class="span6",
                                 checkboxInput(inputId = "normalize",
                                               label = "Normalize read counts",
                                               value = TRUE)),
                             helpText("Scale counts to total aligned reads."),

                             gsub("label class=\"radio\"", "label class=\"radio inline\"",
                                  radioButtons("bias", "Distribution type:",
                                               c("5'" = "left",
                                                 "center" = "middle",
                                                 "3'" = "right"),
                                               "center"))

                             )
                         )

                )
            
            ),
        
        ## Show the caption a line graph of the daily rate and summary of results 
        mainPanel(
            h4(textOutput("caption")),
            tabsetPanel(
                tabPanel("Read Depth",
                         with(tags, 
                              div(class="plot_container",
                                  plotOutput("profile"),
                                  dropButton(inputId = "printmenu1",
                                             label = tags$img(src="navicon.svg"),
                                             choices = c(
                                                 "print"="Print chart",
                                                 "png"="Download PNG image",
                                                 "pdf"="Download PDF document"),
                                             class="pull-right navicon hidden"
                                             ),
                                  uiOutput('viewSlider')
                                  ),
                              ),
                         conditionalPanel(
                             condition="$('div#profile').hasClass('recalculating')",
                             img(src="loading.gif"))
                         ),
                tabPanel("Table", tableOutput("view")),
                tabPanel("NEWS",
                         includeHTML("NEWS.html"))
                )
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
