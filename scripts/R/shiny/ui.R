suppressMessages( library(RMySQL) )
library(shiny)

getOrganisms <- function() {
  query <- paste0('select organism from experiments_tbl group by organism;')

  con = dbConnect(dbDriver("MySQL"), group="morrisdata")
  df = dbGetQuery(con, query)

  return (df$organism)
}


organisms = getOrganisms()
print(organisms)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Miles Per Gallon"),

  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("variable", "Variable:",
                list("Cylinders" = "cyl", 
                     "Transmission" = "am", 
                     "Gears" = "gear")),
    selectInput("org", "Organism:",
                organisms),
    
    checkboxInput("outliers", "Show outliers", FALSE)
  ),

  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    h3(textOutput("caption")),

    plotOutput("mpgPlot")
  )
))
  
