# libraries used. install as necessary

# Time-stamp: <2013-10-21 08:20:40 chris>

  library(shiny)
  library(ggplot2) # graphs
  library(grid)
  library(plyr)  # manipulating data
  library(reshape)
  library(vegan)

## FIXME - this hardcoded path will not be portable.
  source("../shared/R/morrislib.R")
  source("../shared/R/profiling.R")


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
stdPalette <- c("red", "blue")

shinyServer(function(input, output, session) {

    observe({
        message("Inside printmenu observer")
        selection = input$printmenu1
        if (is.null(selection))
            return()
        selection = sub("^([^.]*)-.*", "\\1", selection)
        if (selection == "print") {
            message("	invoke print")
        } else if (selection == "png") {
            message("	invoke png")
        } else if (selection == "pdf") {
            message("	invoke pdf")
        } else if (selection == "stop") {
            message("	invoke stop")
            shiny::stopApp()
        } else {
            message(paste0("	unknown printmenu command - ", selection) )
        }
    })
    
    ## Reactive function to create a list of datasets.  This list may
    ## contain a set of specifications for a database, or it might be a
    ## specification for retrieving data from a file.
    inputDatasets <- reactive({
        if (is.null(input$dataspec))
            return()
        message("in inputDatasets reactive function")
        ## Choose the datasets that meet the user's criteria
        if (input$dataspec == "prostate") {
            control <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+")
            treated <- c("103112_A", "030713_B", "041713_B")
            datasets = c(control=control, treated=treated)
        } else if (input$dataspec == "cll") {
            control = c("061113_A")
            treated = c("061113_B", "061113_C", "061113_D")
            datasets = c(control=control, treated=treated)
        } else if (input$dataspec == "PEO1-RPT-corrUp") {
            control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
            treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "PEO1-RPT-corrUp-coverage.tsv"
        } else if (input$dataspec == "PEO1-RPT-corrDown") {
            control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
            treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "PEO1-RPT-corrDown-coverage.tsv"
        } else if (input$dataspec == "PEO1-RPT-top200") {
            control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
            treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "PEO1-RPT-top200-coverage.tsv"
        } else if (input$dataspec == "SOC") {
            control <- c(
                "Benign_1312_501369.norm", "Benign_1675_224804.norm", "Benign_441_206864.norm",
                "Benign_2638_223151.norm", "Benign_7012_315111.norm", "Benign_7609_361548.norm",
                "Benign_4764_251506.norm", "Benign_1069_501203.norm", "Benign_116_100831.norm",
                "Benign_2186_509428.norm")

            treated <- c(
                "SOC_7637_361542.norm", "SOC_7777_371281.norm", "SOC_12523_494920.norm",
                "SOC_849_206653.norm", "SOC_9547_467919.norm", "SOC_6208_299803.norm",
                "SOC_5991_294171.norm", "SOC_2745_226696.norm", "SOC_5959_278682.norm",
                "SOC_13451_492771.norm")

            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "SOC-merged.tsv"
        } else if (input$dataspec == "test") {
            control = c("C1_RPT.norm", "C2_RPT.norm", "C3_RPT.norm", "C4_RPT.norm")
            treated = c("P1_RPT.norm", "P2_RPT.norm", "P3_RPT.norm", "P4_RPT.norm")
            datasets = c(control = control, treated = treated)
            attr(datasets, "filename") = "test.tsv"
        } else {
            stop(paste0("unknown data specification - ", input$dataspec))
        }

        if (!input$usenorm) {
            attr.saved <- attributes(datasets)
            base <- sub("^(.*)\\..*", "\\1", datasets)
            datasets <- paste0(base, ".depth")
            attributes(datasets) <- attr.saved

            message("datasets:")
            message(paste(datasets, collapse=","))
        }

        return(datasets)
    })


    ## Reactive function to retrieving complete data for a dataset
    dataFromFile <- reactive({
        ## only works for datasets with filename attribute
        datasets = inputDatasets()
        if (is.null(attr(datasets, "filename")))
            return()

        datafile <- attr(datasets, "filename")
        base <- sub("^(.*)\\..*", "\\1",datafile)
        if (file.exists(paste0(base, ".Rda"))) {
            load(paste0(base, ".Rda"))
        } else {
            data <- read.csv(datafile, sep='\t',stringsAsFactors=FALSE)
            data[,'symbol'] <- as.factor(data[,'symbol'])
            data[,'chrm'] <- as.factor(data[,'chrm'])
            save("data", file=paste0(base, ".Rda"))
        }
        # don't resrict what columns are available in the dataset here.
        # do that later when we build the matrix
        message("columns of data:")
        message(paste(names(data), collapse=","))
        return(data)

        ##data <- data[, c("symbol", "transPos", 'exonNumber', datasets)]
    })
        

    # Reactive function for retrieving a list of genes covered in the
    # datasets as specified by inputDatasets.
    geneChoices <- reactive({
        message("in gene choices")
        datasets <- inputDatasets()
        if (is.null(datasets))
            return()

        # if the data is stored in a file, we'll pretty much have to read the whole thing in.
        if (! is.null(attr(datasets, "filename"))) {
            data <- dataFromFile()
            ## extract the list of genes.
            return(levels(data$symbol))
        }

        ## Otherwise the data is stored in a database and we can let
        ## the DB filter out the gene list for us.
        
        ## Get the top expressed genes for the aggregate datasets?
        ## Get all the known genes for the genome?

        ## Try getting the top expressed genes.
        ## Take the datasets, grab the gene list from the dataset
        ## run unique on it.
        ## as.data.frame(table(datasets))
        else if (input$dataspec == "prostate") {
            return(c(c( "Tubb5", "Hist2h2bb", "Nov", "Chgb", "Hist1h1a" ),
                     c("Tgm4", "Pbsn", "Sbp", "Rnase1"),
                     c("Vim", "Itgb1", "Itga1", "Col1a1", "Col1a2")))
        } else  if (input$dataspec == "cll") {
            return(c("RPS2", "RPS23", "EIF3E", "EIF3H", "EIF3I", "EEF1a1", "EIF4EBP2"))
        } else {
            stop(paste0("unknown data specificartion - ", input$dataspec))
        }
        
    })
                
                
    output$geneSelect <- renderUI({
        message("in renderUI of geneselect")
        choices = geneChoices()
        if(is.null(choices))
            return()
        selectInput("geneSelect", "", choices=choices)
    })



    ## Retreive the gene annotation for the currently selected gene in
    ## the currently selected datasets.  To do this we must know the
    ## annotation that was used to assemble the dataset...
    ##
    knowngene <- reactive({
        gene = input$geneSelect
        if(is.null(gene))
            return()

        datasets = inputDatasets()
        if (is.null(datasets))
            return()

        
        # This is a temporary measure useful for datasets retrieved from
        # the database only.  The file-based datasets contain the exon
        # boundaries so no knowngene lookup is necessary.
        #
        # Just return NULL if queried against a file-based dataset
        if (! is.null(attr(datasets, "filename"))) 
            return()
            
        message("in reactive knowngene function")

        ## Otherwise each dataset will tell us which genome it
        ## comes from.  Just use the genome from the first
        ## dataset, as all the datasets should be aligned to the
        ## same genome.
        genome <- morris.getGenome(datasets[[1]])

        kg <- morris.getknowngenes(genome, gene=gene, group=NULL)

        ## FIXME - there may be more than one knowngene (isoform) for
        ## a given common name.  consider what the knowngene data is
        ## being used for and decide if multiple rows in this tables
        ## makes a difference for your purposes.


        return(kg)
    })


    # return the transcript positions of splice junctions in the currently selected gene
    spliceJunctions <- reactive({
        gene = input$geneSelect
        if(is.null(gene))
            return()

        datasets = inputDatasets()
        if (is.null(datasets))
            return()

        message("in reactive spliceJunctions function")

        if (! is.null(attr(datasets, "filename"))) {
            data <- dataFromFile()
            data <- data[data$symbol == gene,]
            isplices <- which(diff(data$exonNumber) != 0)
            splices <- data[isplices, 'transPos']
            return(splices)
        }
        return()
    })
    
    ## Reactive function for retrieving a matrix representing aligned
    ## read positions on a patricular gene for a collection of
    ## datasets.  The datasets might repesent replicated control and
    ## treated conditions, and the gene is one particular gene chosen
    ## from among those in the datasets

    profileMatrix <- reactive({
        gene = input$geneSelect
        if(is.null(gene))
            return()

        datasets = inputDatasets()
        if (is.null(datasets))
            return()

        message("in reactive profileMatrix function")

        if (! is.null(attr(datasets, "filename"))) {
            data <- dataFromFile()
            data <- data[data$symbol==gene,]
            rownames(data) = as.numeric(data$transPos)
            data = subset(data,select=datasets)
            data <- data[order(as.numeric(row.names(data))),]
            mat = as.matrix(t(data))

            if (input$log) {
                # Take the log of each datapoint
                mat <- log10(mat+1)
            }
            
            if (input$centerdata) {
                # center each dataset about its mean
                mat <- t(scale(t(mat), center=TRUE, scale=FALSE))
            }
            
            if (input$renorm) {
                # scale each dataset by the total number of reads in that gene.
                mat <- t(scale(t(mat), center=FALSE, scale=apply(mat,1,sum)))
            }
            
            message("exiting reactive profileMatrix function early")
            return(mat)
        }
        
        genome <- morris.getGenome(datasets[[1]])
        descriptions = morris.fetchinfo(datasets)[,"description", drop=FALSE]
        stats = morris.fetchstats(datasets)

        ## remember the refseq name because that is what identifies each gene in a dataset
        kg = knowngene()
        if (is.null(kg))
            return()
        
        refseq = kg[1,'name']

        ## create an empty matrix to hold the experimental data
        mat <- matrix(0, 0, 50)
        for (dataset in datasets) {
            print(paste0("retrieving data for ", dataset))
            df = morris.getalignments(dataset, refseq)
            attr(df, "dataset") <- descriptions[dataset,"description"]
            ribo.profile = profile(df, kg[refseq,])
            df <- ribo.profile$plotpositions()

            ## count how many reads occur on each position.
            histdata <- hist(df$rposition, breaks=c(1:ribo.profile$transcript()$txLength()), plot=FALSE)
            ## normalize the count at each position by the total RPM
            ## of mapped reads in the dataset.
            scores <- histdata$counts / (stats[dataset,"aligned_count"]/1e6)
            print(paste0("normalization factor for ", dataset, " = ", (stats[dataset,"aligned_count"]/1e6)))
            dim(mat) = c(dim(mat)[1], length(scores))
            mat <- rbind(mat, scores)
            print(paste0("finished data for ", dataset))
        }  ## for each dataset
        rownames(mat) = datasets
        
        message("exiting reactive data function")
        return(mat)
    })


    ## return a distance matrix calcuated from the read profile data.
    ##
    profileDistance <- reactive({
        mat = profileMatrix()
        if (is.null(mat))
            return()
        d <- dist(mat, method=input$distOption)
    })

    # use random permutation of the datasets to calcuate a null distribution
    # and perform an F-test to determine a pseudo p-value.
    ftestData <- reactive({
        ## How likely is the grouping that we see?   If the datasets were partitioned into
        ## groups in a  random fashion, how often would we see a grouping like the observed grouping?

        ## adonis(formula, data, permutations = 999, method = "bray",
        ##        strata = NULL, contr.unordered = "contr.sum",
        ##        contr.ordered = "contr.poly", ...)

        ## There are two ways to call vegan::adonis.  The first used raw data
        ## and the second uses a distance matrix calculated from the raw data.
        ## In the former case, you pass a matrix of raw measurments to
        ## adonis() and it calculates the distance matrix for you using the
        ## algorithm specified by the 'method' argument.  In the second type
        ## of usage, you precompute a distance matrix from the raw data using
        ## any algorithm you like and then pass the distance matrix as the
        ## first argument to adonis().  adonis() can tell which type of usage
        ## you expect by looking at the type of the first argument; if it is
        ## of class 'dist' then you have precomputed a distance matrix,
        ## otherwise you want adonis() to do it for you.
        ##
        ## Why would you want to use adonis() one way versus the other?  It is
        ## usually going to be easier to pass raw data to adonis() and let it
        ## calculate the distance matrix for you.  However, if you want to try
        ## out a distance calculation that is not one of the algorithms
        ## available through adonis() then you will want to precopute your own
        ## distance matrix and pass it in.
        distance <- profileDistance()
        mat = profileMatrix()
        if (is.null(mat))
            return()

        ## get the datasets, as this will tell us which are control
        ## group and which are the experimental group.
        datasets = inputDatasets()
        if (is.null(datasets))
            return()
        treated <- grep("treated", names(datasets))
        control <- datasets[setdiff(1:length(datasets), treated)]
        aVec <- as.factor(sapply(rownames(mat) %in% control, function(x) if (x) {2} else {3}))
        permutations <- 10000

        df <- adonis(distance ~ aVec, permutations=permutations)
        return(df)
    })

    # Create a heading string that is visible at the top of the plot
    # area.  In case youdon't remember why you would do this, this
    # response element is essentially a debugging aid.  It reponds
    # every time the dropdown menu 'printmenu1' is triggered and
    # prints the current time.  Because the time string is always
    # changing, it is easy to see when this has been triggered.
    #
    # To see the output from this response, you must have something
    # like these tags in your ui.R definition:
    #
    # 		h3(textOutput("debug")),
    #
    output$debug <- renderText({
        paste("", Sys.time())
    })
    

    output$coords = renderText({
        gene = input$geneSelect
        if(is.null(gene))
            return()

        ## get the datasets, as this will tell us which are control
        ## group and which are the experimental group.
        datasets = inputDatasets()
        if (is.null(datasets))
            return()
        treated <- grep("treated", names(datasets))
        control <- setdiff(1:length(datasets), treated)

        # get the two dimension points to be drawn on the scatter plot
        fit <- cmdscale(profileDistance(), eig=TRUE, k=2)
        df <- data.frame(fit$points)

        x = (input$click)$x
        y = (input$click)$y
        w <- with(df, which.min( (X1-x)^2 + (X2-y)^2 ))
        sprintf("%s", datasets[[w]])
    })

    ## create 2-D scatter plot of multidimensional sampling data
    output$mdsplot <- renderPlot({
        gene = input$geneSelect
        if(is.null(gene))
            return()

        ## get the datasets, as this will tell us which are control
        ## group and which are the experimental group.
        datasets = inputDatasets()
        if (is.null(datasets))
            return()
        treated <- grep("treated", names(datasets))
        control <- setdiff(1:length(datasets), treated)

        # get the two dismension points to be drawn on the scatter plot
        fit <- cmdscale(profileDistance(), eig=TRUE, k=2)
        df <- data.frame(fit$points)

        ## user can select a standard palette or one that is more
        ## visible to those with R-G color blinkdness.
        palette <- if (input$colorOption) cbPalette else stdPalette

        # make a factor colum that indicates the condition (treated or
        # control) of each dataset
        cond = NULL
        cond[grep("treated", names(datasets))] <- "treated"
        cond[setdiff(1:length(datasets), grep("treated", names(datasets)))] <- "control"
        cond <- factor(cond)

        # add a column indicating what color should be used for each point.
        df <- transform(df, col=palette[cond])
        print((names(datasets)))

        par(mar=c(3, 3, 0.5, 1))  # Trim margin around plot [bottom, left, top, right]

        par(mgp=c(1.5, 0.2, 0))  # Set margin lines; default c(3, 1, 0) [title,labels,line]
        par(xaxs="r", yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

        plot(X2 ~ X1, data=df, col=as.character(col), pch=21, bg=as.character(col),
             xlab="", ylab="", cex=1.5, frame.plot=F, yaxt="n")
        print(par("yaxp"))
        print(paste0(c(min(df$X2), max(df$X2), 6), collapse=", " ))
        ticks = pretty(c(min(df$X2), max(df$X2)), 3)
        print(paste0(ticks, collapse=", "))

        # par(mgp=c(axis.title.position, axis.label.position, axis.line.position))
        mgp <- par("mgp")
        print(paste0(mgp, collapse=", "))
        mgp[2] <- 0.5
        par(mgp=mgp)

        axis(2, at=ticks, labels=T, lwd=0, lwd.ticks=1, lty="solid", las=1, cex.axis=0.9)
        grid()

        # annotate the graph with the pvalue for this gene
        ftest <- ftestData()
        pvalue <- ftest$aov.tab[6][1,1]
        mtext(sprintf("pvalue = %0.4g", pvalue), side=1, line=1, adj=0)

        message("exiting plot")
    })


    ## create line plot for read depth
    output$rdplot <- renderPlot({
        gene <- input$geneSelect
        if(is.null(gene))
            return()

        mat <- profileMatrix()
        if (is.null(mat))
            return()

        datasets <- inputDatasets()
        if (is.null(datasets))
            return()

        message("in renderPlot for read depth plot")

        
        ## get the datasets, as this will tell us which are control
        ## group and which are the experimental group.
        treated = grep("treated", names(datasets))
        control = setdiff(1:length(datasets), treated)

        # This routine will plot a subsection of the data if the
        # current view has been adjusted.  In a normal case, this
        # would just be done with a call to coord_cartesian.  However
        # our data is often so dense that we can get better
        # preformance without sacrificing fidelity by not graphing
        # most of the points.  If the view has been zoomed in then
        # these adjustment need to be made to the points that are
        # still visible.
        #
        # Coord_cartesian is still called so that other annotations,
        # like the splice junctions do not force the view to expand
        # again.
        view = input$viewSlider
        if (is.null(view)) {
            warning("viewSlider returned NULL!")
            view <- c(1,ncol(mat))
        }

        # Bin the data!
        # This is done to make graphing significantly faster without sacrificing much fidelity.
        #
        # If the user has checked the "binning" option, AND if the
        # width of the current plot view is greater than some cutoff
        # (3000 units), then divide the viewable data up into 2000
        # regions and keep the max value within each region.  The X
        # value for each range will be the middle of the range.  When
        # plotting 8 datasets across an 18K bp gene, this technique
        # reduced the number of points from (18K x 8) to (2K x 8).
        #
        # How were the magic numbers 2000 and 3000 chosen?  2000 was
        # chosen because the graph still looked good at that
        # resolution on my macbook air screen.  3000 because I didn't
        # see much point in binning if the complete data was only a
        # small fraction wider than the target resolution.
        x <- seq(view[[1]], view[[2]])
        if (input$bindata && length(x) > 3000) {
            bpt <-  pretty(x, n=2000)
            x.cut <- cut(x, bpt)
            mat <- t(apply(mat[,x], 1, function(r) sapply(split(r, x.cut), max)))
            x.mean <- sapply(split(x, x.cut), mean)
            colnames(mat) <- x.mean
        } else {
            mat <- mat[,x]
            colnames(mat) <- x
        }
        
        ## user can select a standard palette or one that is more
        ## visible to those with R-G color blindness.
        colorPalette <- if (input$colorOption) cbPalette else stdPalette
        yaxis_label  <- if (input$log) "read depth log2(normalized to RPM)" else "read depth (normalized to RPM)"
            

        gg <- ggplot(melt(mat), aes(name="", x=X2,
                                    y=value,
                                    group=X1), environment = environment())
        gg <- gg + theme_bw()
        gg <- gg + theme(legend.title = element_text(size = 16, face = "bold"),
                         legend.text = element_text(size = 14, face = "bold"),
                         legend.position="top",
                         legend.direction="horizontal")
        gg <- gg + theme(legend.key = element_rect(size = 0.5, linetype="blank"))
        gg <- gg + theme(panel.border = element_blank())
        gg <- gg + theme(plot.margin = unit(c(0,0,0,0), "cm"))

        gg <- gg + ylab("read depth (normalized to RPM)")
        gg <- gg + xlab("Transcript position")
        gg <- gg + ylim(-max(mat),max(mat))
        gg <- gg + scale_colour_manual(name=gene, values=colorPalette,
                                       labels=c("Control", "Treated"))
        gg <- gg + geom_line(data=melt(mat[treated,,drop=FALSE]), aes(color=colorPalette[[1]]))
        gg <- gg + geom_line(data=melt(-mat[control,,drop=FALSE]), aes(color=colorPalette[[2]]))

        gg <- gg + coord_cartesian(xlim = view)

        if (input$showSplices) {
            splices <- spliceJunctions()
            if (!is.null(splices)) {
                gg <- gg + geom_vline(xintercept = splices, alpha=.25,  color="red", linetype="dashed")
            }
        }
        print(gg)
    })

    
    output$viewSlider <- renderUI({
        mat <- profileMatrix()
        if (is.null(mat))
            return()
        
        sliderInput(inputId = "viewSlider",
                    label=" ",
                    min = 1, max = ncol(mat), step = 1,
                    value = c(1,ncol(mat)))
    })

    
    ## create summary data for each subject 
    output$view <- renderTable({
        mat = data()
        if (is.null(mat))
            return()
        mat
    })
    
    ## make data downloadable
    output$downloadData <- downloadHandler(
        ## filename = function() { paste(input$data, '.csv', sep='') }, in tutorial to distinguish files trickier with my work
        filename = function() { paste('results.csv', sep='') },
        content = function(file) {
            write.csv(data(), file)
        }
        )
    
})
