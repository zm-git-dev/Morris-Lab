require(shiny)
source("~/Morris-Lab/scripts/R/morrislib.R")
source("~/Morris-Lab/scripts/R/profiling/profiling.R")

require(Biostrings)
require(genomes)



Sys.setenv(email='cswarth@gmail.com')


nuc2aa <- function(gene) {

    ## use bioconductor - ShortRead package to read the fasta format.
    ## http://stackoverflow.com/a/9331406/1135316
    ## http://www.bioconductor.org/packages/release/bioc/html/ShortRead.html

    ## looks like BioStrings is a better package to read fasta strings.
    ## use the genomes package the fetch data from NCBI
    ## http://www.bioconductor.org/packages//2.10/bioc/manuals/genomes/man/genomes.pdf
    ## Results must be saved to a temp file
    ## http://www.biostars.org/p/75700/#75709
    tmp = tempfile()
    efetch(gene, showURL=FALSE, db="nucleotide", retmode="text", rettype="fasta", destfile=tmp, complexity=3, strand=1)
    records = readBStringSet(tmp)
    mrna = as(records[[1]], "DNAString")
    protein = as(records[[2]], "AAString")

    ## http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter7.html
    ## Find all occurrences of codon "atg" in sequence "mrna"
    match.atg <- matchPattern("atg", mrna)

    ## extract the subsequence starting at the first start codon.
    s2b <- subseq(subject(match.atg), start=start(match.atg)[1], end=length(subject(match.atg)))

    ## translate to codons and to amino acid sequence.
    codons = as.character(codons(s2b))
    aa = strsplit(as.character(translate(s2b)), '')[[1]]

    ## extract the codons and amino acid sequence up to the first stop codon.
    n = match('*', aa)
    data.frame(codon=codons[1:n], aa=aa[1:n])

}

# Define server logic for random distribution application
shinyServer(function(input, output) {

  # Reactive expression to generate the requested alignments. This is 
  # called whenever the inputs change. Some of the renderers defined 
  # below then  use the value computed from this expression
  alignments <- reactive({
      gene <- input$gene     # either common name or refseq name
      dataset <- input$dataset
      genome = morris.getGenome(dataset)

      kg <- morris.getknowngenes(genome, gene=gene, group=NULL)
      rownames(kg) <- kg$name
      ## should check that we get a unique gene!
      df = morris.getalignments(dataset, rownames(kg)[1])

      profile(df,  kg[1,])
  })

  codons <- reactive({
      nuc2aa(input$gene)
  })
  
  # Generate a plot of the data. Also uses the inputs to build the 
  # plot label. Note that the dependencies on both the inputs and
  # the 'data' reactive expression are both tracked, and all expressions 
  # are called in the sequence implied by the dependency graph
  output$plot <- renderPlot({
      plot(alignments(), minlen=28, units=input$units)
  })

  # Generate a summary of the data
  output$summary <- renderPrint({
    summary(alignments())
  })

  # Generate an HTML table view of the data
  ##output$table <- renderTable({
  output$table <- renderTable({
      codons()
  })
})
