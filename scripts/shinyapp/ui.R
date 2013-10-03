library(shiny)
source("~/Morris-Lab/scripts/R/morrislib.R")

datasets = morris.datasets(organism="Mouse")
df.info = morris.fetchinfo(datasets)
print(df.info)
choices=rownames(df.info)
names(choices)=paste0(df.info$experiment, " ", df.info$tissue, " ", df.info$description)

# Define UI for random distribution application 
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Profiles"),

  ## Sidebar with controls to select the random distribution type
  ## and number of observations to generate. Note the use of the br()
  ## element to introduce extra vertical spacing

  sidebarPanel(
               selectInput("dataset", "Choose a dataset:", 
                           choices = choices),
               br(),
               
               textInput("gene", "Gene:", value = "NM_007409"),

    radioButtons("units", "Distribution type:",
                 list("Nucleotides" = "nucleotides",
                      "Codons" = "aa")),
    br(),

    sliderInput("n", 
                "Number of observations:", 
                 value = 500,
                 min = 1, 
                 max = 1000)
  ),

  # Show a tabset that includes a plot, summary, and table view
  # of the generated distribution
  mainPanel(
    tabsetPanel(
      tabPanel("Plot", plotOutput("plot")), 
      tabPanel("Summary", verbatimTextOutput("summary")), 
      tabPanel("Table", tableOutput("table"))
    )
  )
))
