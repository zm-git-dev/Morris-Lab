
library(grid)


morris.getalignments <- function(dataset, gene, group=NULL) {
    result <- tryCatch({
        drv <- dbDriver("MySQL")
        if (is.null(group)) {
            con <- dbConnect(drv, group="remote")
            # con <- dbConnect(drv, user="readonly", password="readonly", dbname="morris", host="localhost")
        } else  {
            con <- dbConnect(drv, group=group)
        }
        if (is.null(con)) {
            stop(paste0("Could not connect to database: ", e$message));
        }

        ## SQL query for retrieving the per-gene alignment count of a single dataset.
        ## These individual tables will be joined together based on the gene name
        ## as a common key.   Oddly eough it is faster to do this in R than it is
        ## in SQL.
        query <- paste0("select feature,position,length from ",
                        "new_alignments_tbl a join datasets_tbl d on a.dataset_id=d.id ",
                        "where d.name like '", dataset, "%' and a.feature in (",
                        paste0("'",gene,"'", collapse=","), ")")

        df <- dbGetQuery(con, query)
        if (nrow(df) == 0) {
            print(paste0("Error during database query: no data in dataset '", dataset,"'"))
            print("Available datasets:")
            print(morris.datasets(group));
            stop("No Data");
        }

        query <- paste0("select genome as '", dataset,"' from datasets_tbl d ",
                        "where d.name like '", dataset, "%'")
        gf <- dbGetQuery(con, query)
            if (nrow(df) == 0) {
                print(paste0("Error during database query: no dataset named '", dataset,"'"))
                print("Available datasets:")
                print(morris.datasets(group))
                stop("No Data")
            }
        attr(df, "genome") <- gf[1,1]
        attr(df, "dataset") <- dataset
    }, finally = {
        if (exists("con")) 
          dbDisconnect(con)
        ##if (exists("drv")) 
        ##  dbUnloadDriver(drv)
    })

    return(df)
}

## This awesome flavor of OOP in R using closures is taken from a
## stackoverflow posting I stumbled across:
## http://stackoverflow.com/a/15245568/1135316
## Add a few lines and it is now an S3 style class.
transcript = function(gdata) {
    starts <- as.numeric(strsplit(gdata$exonStarts,",")[[1]])
    ends <- as.numeric(strsplit(gdata$exonEnds,",")[[1]])
    elen <- ends - starts
    transcriptLength <- sum(elen)

    length = function() transcriptLength
    cdsLength = function() cdsEnd()-cdsStart()
    cdsStart = function() rpos(gdata$cdsStart)
    cdsEnd = function() rpos(gdata$cdsEnd)
    txStart = function() rpos(gdata$txStart)
    txEnd = function() rpos(gdata$txEnd)
    isCoding = function() (gdata$cdsEnd != gdata$cdsStart)
    name = function() (gdata$name)
    name2 = function() (gdata$name2)
    rpos = function(pos) {
        ## Given a position on a chromosome, the position within transcript is
        ## equal to the sum of all exons that end before the
        ## position, PLUS the beginning portion of the exon that contains the position.

        positionInTranscript = 0
        ## there are two primary situations:
        ## 1) the position is to the right of the start of the transcript
        ## 2) the position is to the left of the start of the transcript
        if (pos >= starts[1]) {
            ## positioned to the right of the beginning of the gene
            if (any(ends<pos)) {
                positionInTranscript = sum( elen[ends<pos] )
            }
            if (any(ends > pos)) {
                positionInTranscript = positionInTranscript + (pos - starts[ends > pos][1] + 1)
            } else {
                positionInTranscript = positionInTranscript + (pos - tail(ends,1))
            }
        }
        
        if (gdata$strand == '-') {
            positionInTranscript = transcriptLength - positionInTranscript
        }
        return(positionInTranscript)
    }
    
    exported = list(
      length=length,
      cdsLength=cdsLength,
      cdsStart=cdsStart,
      cdsEnd=cdsEnd,
      txStart=txStart,
      txEnd=txEnd,
      isCoding=isCoding,
      name=name,
      name2=name2,
      rpos=rpos
      )
    class(exported) <- "transcript"
    invisible(exported)
}

print.transcript <- function(object) {
    cat("transcript ", object$name2(), ", ", object$length(), "\n")
}

plot.transcript <- function(self, xlim=NULL) {
    if (is.null(xlim))
      xlim <- c(1, self$length())
    clip(xlim[1], xlim[2], 0, 100)

    cdsHeight = 5   # height of coding region, expressed as percentage of plot height
    segments(xlim[1], 100-cdsHeight,  xlim[2],  100-cdsHeight, lwd=5, lend="butt")
    if (self$isCoding()) {
        ## Label the endpoints of the coding region and place the gene
        ## name in the middle of the coding region.

        rect(self$cdsStart(), 100-2*cdsHeight, self$cdsEnd(), 100, col='blue')
        text((self$cdsStart() + self$cdsEnd())/2, 100-2*cdsHeight-5,
             adj=c(.5,0.5), self$name2())
        text(self$cdsStart(), 100-2*cdsHeight-5, adj=c(.5,0.5), as.character(self$cdsStart()))
        text(self$cdsEnd(), 100-2*cdsHeight-5, adj=c(.5,0.5), as.character(self$cdsEnd()))
    } else {
        ## There is no coding region - this is a non-coding gene.
        ## Label the endpoints of the transcript and label the gene in the middle.
        text( self$length()/2, 100-2*cdsHeight-5,
             adj=c(.5,0.5), self$name2())
        text(1, 100-2*cdsHeight-5, adj=c(.5,0.5), as.character(0))
        text(self$length(), 100-2*cdsHeight-5, adj=c(.5,0.5), as.character( self$length() ))
    }
}

profile <- function(df, transcript, xlim=NULL) {
    usr = par("usr")
    
    ## Calculate the relative position of each alignment w.r.t. start of transcript.
    positions = sapply(X=df$position, FUN=transcript$rpos)

    ## save the graphical environent.
    plt <- par("plt")

    ## prepare the graph the histogram in the top half of the figure.
    par(plt=c(plt[1], plt[2], .5, 1))

    ## Calculate the histogram bins and maximum count across all bins.
    ##
    if (is.null(xlim))
      xlim <- c(1, transcript$length())

    x <- positions[!((positions < xlim[1]) | (positions > xlim[2]))]
    binsize <- 1   ## for small transcripts with few reads, use a bin size of one.
    if (length(x) > 50) 
        binsize <- 3  ## otherwise use a bin size of three
    histdata <- hist(x, breaks=seq(1,transcript$length(), binsize), plot=F)
    maxy <- max(histdata$count)

    plot(histdata$mids[histdata$counts != 0],histdata$counts[histdata$counts != 0], type='h',
         xlim=xlim, ylim=c(0,1.25*maxy), 
         lwd=3, lend=2,ylab='',xlab='',main='',xaxt='n',bty='n')

    # The usr coordinates appear to change after plot(), so fetch them again.
    tmp = par("usr")
    text(xlim[1], tmp[4], attr(df, "dataset"), adj=c(0,1.2), new=TRUE)
    text(xlim[2], tmp[4], paste0(length(x), " total reads"), adj = c( 1, 1.2 ), new=TRUE)

    ## start a new plot without clearing the current one.
    par(new=TRUE)
    plot.new()

    ## draw the transcript in the lower half of the plot area
    par(plt=c(plt[1], plt[2], 0, .5))
    par(usr=c(tmp[1], tmp[2], 0, 100))
    plot(transcript, xlim=xlim)
    par(usr=usr)

}


.pardefault <- par(no.readonly = T)

#gene="NR_029560"   # Mir150
gene='NM_001005419'	# 'Ado'  single exon reverse
gene='NM_013541'	# 'Gstp1'
gene='NM_016978'	# 'Oat'
gene='NM_011434'	# 'Sod1'
gene='NM_144903'	# Aldob
gene='NM_011044'	# 'Pck1'
gene="NM_007409"   # ADH1
#gene='NM_TEST'		# 'test'

plot.new()

    df = morris.getalignments("113010_A", gene)
    ##df=data.frame(position=c(2,2,2,3,4,4,5,5,5,5,7,7,9,9,9), length=22)

    kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
    ##kg <- data.frame(genome="mm9", name=c("gene1","gene2"), strand='+', txStart=1, txEnd=20, cdsStart=c(5,6), cdsEnd=c(8,9), exonCount=1, exonStarts="1,", exonEnds="20,", name2=c("test","test2"),stringsAsFactors = F)
    rownames(kg) <- kg$name

    gobj = transcript(kg[gene,])

    profile(df, gobj)
##    profile(df, gobj, xlim=c(80,200))

par(.pardefault)
