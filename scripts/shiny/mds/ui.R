

require(shiny)

# Create the Shiny binding object for our component, and register it:
#
dropButton <- function(inputId, 
                        label, 
                        choices, 
                        selected = NULL, 
                        multiple = FALSE) {

    # resolve names
    choices <- shiny:::choicesWithNames(choices)
    
    menuList <- tags$ul(class = "dropButton dropdown-menu pull-right", id=inputId)
    
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
                 tags$a(label, class="btn btn-small dropdown-toggle",
                        "data-toggle"="dropdown", href="#"),
                 menuList))

    print(dropTag)
    dropTag
}


shinyUI(
    pageWithSidebar(
        
        ## Application title
        headerPanel(list("Characterizing ribosome profiles", img(src="FHlogo.gif", width=200, class="pull-right")),
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
            ## h3(textOutput("caption")),
            tags$head( tags$link(rel="stylesheet", type="text/css", href="app.css")),
            tags$head( tags$script(src="app.js")),
            tabsetPanel(
                tabPanel("Read Depth",
                         ## bootstrap dropdown menu code shamelessly stolen from
                         ## http://stackoverflow.com/a/13998987/1135316
                         ##
                         with(tags, 
                              div(class="plot_container", plotOutput("rdplot"),
                                  dropButton(inputId = "printmenu1",
                                             label = tags$img(src="navicon.svg"),
                                             choices = c(
                                                 "print"="Print chart",
                                                 "png"="Download PNG image",
                                                 "pdf"="Download PDF document")
                                             ),
                                  ## div(class="dropdown btn-group navicon",
                                  ##     a(img(src="navicon.svg"),
                                  ##       class="btn btn-small dropdown-toggle",
                                  ##       "data-toggle"="dropdown", href="#"),
                                  ##     ul(class="dropdown-menu pull-right", id="printmenu1",
                                  ##        li(a("Print chart", href="#")),
                                  ##        li(class="divider"),
                                  ##        li(a("Download PNG image", id="alertMe")),
                                  ##        li(a("Download PDF Document", href="#"))
                                  ##        )),
                                  span()
                                  )
                              )
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
                                         li(a("Print chart", href="#")),
                                         li(class="divider"),
                                         li(a("Download PNG image", href="#")),
                                         li(a("Download JPEG image", href="#")),
                                         li(a("Download PDF Document", href="#"))
                                         ))))
                         ), 
                tabPanel("Table", tableOutput("view"))
                )
            ) ## end-mainpanel
        )  ## end-pagewithsidebar
    )  ## end-shinyui
