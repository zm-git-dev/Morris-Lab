 # libraries used. install as necessary

  library(shiny)
  library(RJSONIO) # acquiring and parsing data
  library(ggplot2) # graphs
  library(plyr)  # manipulating data
  library(lubridate) #dates
  library(stringr)

## FIXME - this hardcoded path will not be portable.
  source("~/Morris-Lab/scripts/R/morrislib.R")
  source("~/Morris-Lab/scripts/R/profiling/profiling.R")

trim.leading <- function (x)  sub("^\\s+", "", x)

shinyServer(function(input, output) {
  
  
  data <- reactive(function() {  


      ## Choose the datasets that meet the user's criteria
      if (input$dataspec == "prostate") {
          control <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+")
          treated <- c("103112_A", "030713_B", "041713_B")
      } else if (input$dataspec == "cll") {
          control = c("061113_A")
          treated = c("061113_B", "061113_C", "061113_D")
      } else {
          stop(paste0("unknown data specificartion - ", input$dataspec))
      }
      datasets <- c(treated,control)

      genome <- morris.getGenome(datasets[[1]])
      descriptions = morris.fetchinfo(datasets)[,"description", drop=FALSE]
      stats = morris.fetchstats(datasets)

      gene = input$gene
      kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
      rownames(kg) <- kg$name
      stopifnot(nrow(kg) == 1)

      ## remember the refseq name because that is what identifies each gene in a dataset
      refseq = kg[1,'name']

      ## create an empty matrix to hold the experimental data
      mat <- matrix(0, 0, 50)
      for (dataset in datasets) {
          df = morris.getalignments(dataset, refseq)
          attr(df, "dataset") <- descriptions[dataset,"description"]
          ribo.profile = profile(df, kg[refseq,])
          df <- ribo.profile$plotpositions()
          xlim <- c(1, ribo.profile$transcript()$txLength())

          ## count how many reads occur on each position.
          histdata <- hist(df$rposition, breaks=c(1:ribo.profile$transcript()$txLength()), plot=FALSE)
          scores <- histdata$counts * (mean(stats[datasets, "raw_count"])/stats[dataset,"raw_count"])
          scores <- histdata$counts
          dim(mat) = c(dim(mat)[1], length(scores))
          mat <- rbind(mat, scores)
      }  ## for each dataset
      rownames(mat) = datasets
      
   return(mat)
  })
  
  
  # Create a heading based on range of dates selected for printing as a caption
  output$caption <- renderText(function() {
      paste(input$dataspec, ":", input$gene)
  })
  
 
  # create plot for linear and log scales
  output$plot <- renderPlot(function() {
      mat = data()
      d <- dist(mat)
      fit <- cmdscale(d, eig=TRUE, k=2)
      x <- fit$points[,1]
      y <- fit$points[,2]
      gg <- ggplot(data.frame(x=x, y=y), aes(name="", x=x, y=y, label=rownames(mat)),
                   environment = environment())
      gg <- gg + theme_bw()
      gg <- gg + theme(legend.position="top",legend.title=element_blank(), legend.text = element_text(colour="blue", size = 14, face = "bold"))
      gg <- gg + geom_point() 
      gg <- gg + geom_text()
      gg <- gg + scale_x_continuous(expand = c(.2,0))
      print(gg)

    })
  
 
  # create summary data for each subject 
  output$view <- renderTable(function() {
    summary(data())
  })
  
  # make data downloadable

  output$downloadData <- downloadHandler(
   # filename = function() { paste(input$data, '.csv', sep='') }, in tutorial to distinguish files trickier with my work
    filename = function() { paste('results.csv', sep='') },
    content = function(file) {
      write.csv(data(), file)
    }
  )
  
})
