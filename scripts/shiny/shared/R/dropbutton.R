
require(shiny)

## Code to hook up a bootstrap dropdown button to shiny.

## Copyright (C) 2013 by Chris Warth

timestamp <- "dropbutton.R 2013/09/12 11:34:02 chris"



# Create the Shiny binding object for our component, and register it:
#
dropButton <- function(inputId, 
                       label, 
                       choices, 
                       selected = NULL, 
                       multiple = FALSE
                       ) {
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
        singleton(tags$head(tags$script(src = "dropbutton.js"))),
        tags$div(class = "dropdown btn-group navicon",
                 type="navicon",
                 tags$a(label, class="btn btn-small dropdown-toggle",
                        "data-toggle"="dropdown", href="#"),
                 menuList))

    message(dropTag)
    dropTag
}
message("inside dropbutton.R")

## Local Variables:
## time-stamp-pattern: "10/timestamp <- \"%f %:y/%02m/%02d %02H:%02M:%02S %u\""
## End: 
