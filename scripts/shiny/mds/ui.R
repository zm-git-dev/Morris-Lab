

library(shiny)

# Create the Shiny binding object for our component, and register it:
#
dropButton <- function(inputId, 
                        label, 
                        choices, 
                        selected = NULL, 
                        multiple = FALSE) {

    # resolve names
    choices <- shiny:::choicesWithNames(choices)
    
    menuList <- tags$ul(class = "dropdown-menu pull-right", id=inputId)
    
    # Create tags for each of the options
    ids <- paste0(inputId, seq_along(choices))
    liId <- 1
    for (choice in names(choices)) {
        thisId <- paste("menu", inputId, liId, sep="-")
        liId <- liId + 1
        
        liTag <- tags$li(tags$a(choices[choice],
                                id=thisId, href="#"))
        menuList <- tagAppendChild(menuList, liTag)
    }

    dropTag <- tagList(
        singleton(tags$head(tags$script(src = "app.js"))),
        tags$div(class = "dropdown btn-group navicon",
                 type="navicon",
                 tags$a(label,
                        class="btn btn-small dropdown-toggle",
                        "data-toggle"="dropdown", href="#"),
                 menuList))

    print(dropTag)
    tagList(shiny:::controlLabel(inputId, label), dropTag)
}


shinyUI(
    pageWithSidebar(
        
        ## Application title
        headerPanel("Characterizing ribosome profiles"),
        
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


                         div(class="span6", checkboxInput(inputId = "log",
                                 label = "log10 scale", value = FALSE)),
                         
                         helpText("Use log scale if compared searches are significantly different"),
                         downloadButton('downloadData', 'Download Output as csv')

                         ),
                tabPanel("Options",
                         helpText((paste0(
                             "Select options to control the analysis of ribosome profiles."))),
                         wellPanel(
                             selectInput(inputId = "distOption",
                                         label = "Distance metric:",
                                         choices = c(
                                             "euclidean", "maximum", "manhattan", "canberra",
                                             "binary", "minkowski")
                                         ),
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
            h3(textOutput("caption")),
            tags$head( tags$link(rel="stylesheet", type="text/css", href="app.css")),
            tags$head( tags$script(src="app.js")),
            tabsetPanel(
                tabPanel("Read Depth",
                         dropButton(inputId = "testing",
                                    label = tags$img(src="navicon.svg"),
                                    choices = c(
                                        "Stacked reads"="stacked",
                                        "Centered reads"="centered")
                                    ),
                         ## bootstrap dropdown menu code shamelessly stolen from
                         ## http://stackoverflow.com/a/13998987/1135316
                         ##
                         with(tags, 
                              div(class="plot_container", plotOutput("rdplot"),
                                  div(class="dropdown btn-group navicon",
                                      a(img(src="navicon.svg"), class="btn btn-small dropdown-toggle", "data-toggle"="dropdown", href="#"),
                                      ul(class="dropdown-menu pull-right", id="printmenu1",
                                         li(a("Print chart", href="#", target="_blank")),
                                         li(class="divider"),
                                         li(a("Download PNG image", id="alertMe", target="_blank")),
                                         li(a("Download PDF Document", href="#" target="_blank"))
                                         ))))
                         ),
                
                
                tabPanel("MDS", 
                         ## bootstrap dropdown menu code shamelessly stolen from
                         ## http://stackoverflow.com/a/13998987/1135316
                         ##

                         with(tags, 
                              div(class="plot_container", plotOutput("mdsplot"),
                                  div(class="dropdown btn-group navicon",
                                      a(img(src="navicon.svg"), class="btn btn-small dropdown-toggle", "data-toggle"="dropdown", href="#"),
                                      ul(class="dropdown-menu pull-right", id="printmenu2",
                                         li(a("Print chart", href="#", target="_blank")),
                                         li(class="divider"),
                                         li(a("Download PNG image", href="#", target="_blank")),
                                         li(a("Download JPEG image", href="#", target="_blank")),
                                         li(a("Download PDF Document", href="#", target="_blank"))
                                         ))))
                         ), 
                tabPanel("Table", tableOutput("view"))
                )
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
