
## Code to hook up a bootstrap dropdown button to shiny.

## Copyright (C) 2013 by Chris Warth

timestamp <- "dropbutton.R 2013/09/24 12:06:58 chris"

require(shiny)

## options for dropbuttons
.dropButton.js = "js/dropbutton.js"  # where to read dropbutton javascript from 

##' Create a Bootstrap drop button control.  A drop button responds to
##' mouse clicks (or other navigation?) with a menu from which items
##' may be selected.
##'
##' Be aware that the drop downbutton from v2.1.0 of Bootstrap is
##' broken; it is not automaically cleared after a selection is made.
##' V2.3.2 of Bootstrap does not show the behavior and the menu
##' disappears as one would expect.  I have tried Shiny on top of
##' v2.3.2 and it seems to function OK.  The trick is to fool Shiny
##' into using that version of Bootstrap.
##' @param inputId 
##' @param label 
##' @param choices 
##' @param selected 
##' @param multiple 
##' @return 
##' @author Chris Warth
##' @export
dropButton <- function(inputId, 
                       label, 
                       choices, 
                       selected = NULL, 
                       multiple = FALSE,
                       class=""
                       ) {
    # resolve names
    choices <- shiny:::choicesWithNames(choices)
    
    menuList <- tags$ul(class = "dropButton dropdown-menu", id=inputId)
    
    # Create tags for each of the options
    ids <- paste0(inputId, seq_along(choices))
    liId <- 1
    for (choice in names(choices)) {
        thisId <- paste("menu", inputId, liId, sep="-")
        liId <- liId + 1
        
        liTag <- tags$li(tags$a(choices[choice],
                             id=thisId, "data-value"=choice, href="#", class="shiny-download-link"))
        menuList <- tagAppendChild(menuList, liTag)
    }

    dropTag <- tagList(
        singleton(tags$head(tags$script(src = .dropButton.js))),
        tags$div(class = paste("dropdown btn-group",class),
                 tags$a(label, class="btn btn-small dropdown-toggle",
                        "data-toggle"="dropdown", href="#"),
                 menuList))
}

## Local Variables:
## time-stamp-pattern: "10/timestamp <- \"%f %:y/%02m/%02d %02H:%02M:%02S %u\""
## End: 
