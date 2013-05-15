
##library(grid)


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

##
## I've since realize that this is just the S3 style of classes in R.
## Add a few lines and it is now an S3 style class.
transcript = function(gdata) {
    starts <- as.numeric(strsplit(gdata$exonStarts,",")[[1]])
    ends <- as.numeric(strsplit(gdata$exonEnds,",")[[1]])
    elen <- ends - starts
    transcriptLength <- sum(elen)
    ##print(transcriptLength)

    cdsStart = function() rpos(if (gdata$strand == '+') { gdata$cdsStart } else { gdata$cdsEnd })
    cdsEnd = function() rpos(if (gdata$strand == '+') { gdata$cdsEnd-1 } else { gdata$cdsStart+1 })
    cdsLength = function() cdsEnd()-cdsStart()+1
    peptideLength = function() cdsLength()/3
    txStart = function()  rpos(if (gdata$strand == '+') { gdata$txStart } else { gdata$txEnd })
    txEnd = function() rpos(if (gdata$strand == '+') { gdata$txEnd-1 } else { gdata$txStart+1 })
    length = function() txEnd()-txStart()+1
    isCoding = function() (gdata$cdsEnd != gdata$cdsStart)
    name = function() (gdata$name)
    name2 = function() (gdata$name2)
    rpos = function(pos) {
        ## Given a position on a chromosome, the position within transcript is
        ## equal to the sum of all exons that end before the
        ## position, PLUS the beginning portion of the exon that contains the position.

        offset = 0
        ## there are two primary situations:
        ## 1) the position is to the right of the start of the transcript
        ## 2) the position is to the left of the start of the transcript
        if (pos >= starts[1]) {
            ## positioned to the right of the beginning of the gene
            if (any(ends<=pos)) {
                offset = sum( elen[ends<=pos] )
            }
            if (any(ends > pos)) {
                offset = offset +  (pos - starts[ends > pos][1] + 1)
            } else {
                offset = offset + (pos - tail(ends,1)+1)
            }
        } else {
            offset = pos-starts[1]
        }
        
        if (gdata$strand == '-') {
             offset = transcriptLength - offset+2
        }
        return(offset)
    }
    
    stopifnot(txStart() == 1)
    
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
      rpos=rpos,
      peptideLength=peptideLength
      )
    class(exported) <- "transcript"
    invisible(exported)
}

print.transcript <- function(object) {
    cat("transcript ", object$name2(), ", ", object$length(), "\n")
}

plot.transcript <- function(self, xlim=NULL, units="nucleotide") {
    if (units == "aa") {
        if (is.null(xlim)) 
          xlim <- c(self$cdsStart(), self$cdsEnd())
    } else {
        if (is.null(xlim))
          xlim <- c(1, self$length())
    }
    
    cdsHeight = 10   # height of coding region, expressed as percentage of plot height
    if (self$isCoding()) {
        ## Label the endpoints of the coding region and place the gene
        ## name in the middle of the coding region.
        if (units == "aa") {
            cdslim = c(max(xlim[1], 1), min(xlim[2], self$peptideLength()))
        } else {
            cdslim = c(max(xlim[1], self$cdsStart()), min(xlim[2], self$cdsEnd()))
        }
        text((cdslim[1] + cdslim[2])/2, 100-2*cdsHeight-5,
             adj=c(.5,1.0), self$name2())
        text(cdslim[1], 100-2*cdsHeight-5, adj=c(.5,1.0), as.character(cdslim[1]))
        text(cdslim[2], 100-2*cdsHeight-5, adj=c(.5,1.0), as.character(cdslim[2]))
        clip(xlim[1], xlim[2], 0, 100)
        segments(xlim[1], 100-cdsHeight,  xlim[2],  100-cdsHeight, lwd=5, lend="butt")
        rect(cdslim[1], 100-2*cdsHeight, cdslim[2], 100, col='blue')
    } else {
        ## There is no coding region - this is a non-coding gene.
        ## Label the endpoints of the transcript and label the gene in the middle.
        text((xlim[1]+xlim[2])/2, 100-2*cdsHeight-5,
             adj=c(.5,1.0), self$name2())
        text(xlim[1], 100-2*cdsHeight-5, adj=c(.5,1.0), as.character(xlim[1]))
        text(xlim[2], 100-2*cdsHeight-5, adj=c(.5,1.0), as.character(xlim[2]))
        clip(xlim[1], xlim[2], 0, 100)
        segments(xlim[1], 100-cdsHeight,  xlim[2],  100-cdsHeight, lwd=5, lend="butt")
    }
}

profile <- function(df, gene.data) {
    gene <- transcript(gene.data)
    transcript <- function() gene
    alignments <- function() df
    exported = list(transcript=transcript, alignments=alignments )
    class(exported) <- "profile"
    invisible(exported)
}

plot.profile <- function(self, xlim=NULL, units="nucleotide", bias="middle", minlen=NULL) {
    usr = par()$usr
    plt = par()$plt

    df <- self$alignments()
    transcript <- self$transcript()
    
    ## Calculate the relative position of each alignment w.r.t. start of transcript.
    ## position may be calculated with respect to the middle of the read or
    ## the 5'-end
    if (bias == "middle") {
        x = sapply(X=df$position+df$length/2.0, FUN=transcript$rpos)
    } else {
        x = sapply(X=df$position, FUN=transcript$rpos)
    }
    ## if the user specified a minimum read length, discard shorter reads.
    if (!is.null(minlen))
        x = x[df$len >= minlen]

    ## if the user is zooming in on a portion of the transcript, set the
    ## limits of the horizontal axis accordingly.   Otherwise the limits
    ## are determined by the length of the transcript.
    if (units == "aa") {
        if (is.null(xlim)) {
            xlim = c(transcript$cdsStart(), transcript$cdsEnd())
        } 
        ## convert nucleotide positions to codon position
        x = (x - transcript$cdsStart())/3
    } else {
        if (is.null(xlim))
          xlim <- c(1, transcript$length())
    }

    
    ## discard any alignments that fall outside the horizontal limits.
    ## ( why?  won't the plot just truncate any data outside the range? )
    is.between <- function(x,lim) {
        (x > lim[1]) & (x < lim[2])
    }
    x <- x[is.between(x,xlim)]
    
    ## Draw the histogram in the top 2/3 of the plot area.
    par(plt=c(plt[1], plt[2], plt[3]+(plt[4]-plt[3])/3, plt[4]))

    ## Calculate the histogram bins and maximum count across all bins.
    ##
    if (units == "aa") {
        ## for small transcripts with few reads, use a bin size of one.
        breaks=seq(1, transcript$peptideLength(), 1)
    } else {
        ## otherwise use a bin size of three
        breaks=seq(1,transcript$length(), 3)
    }
    histdata <- hist(x, breaks=breaks, plot=F)
    maxy <- max(histdata$count)

    plot(histdata$mids[histdata$counts != 0],histdata$counts[histdata$counts != 0], type='h',
         xlim=xlim, ylim=c(0,1.25*maxy), 
         lwd=3, lend=2,ylab='',xlab='',main='',xaxt='n',bty='n')

    # The usr coordinates appear to change after plot(), so fetch them again.
    tmp = par("usr")
    text(xlim[1], tmp[4], attr(df, "dataset"), adj=c(0,1.2), new=TRUE)
    text(xlim[2], tmp[4], paste0(length(x), " total reads"), adj = c( 1, 1.2 ), new=TRUE)

    ## start a new plot without clearing the current one.
    ## I don't understand why this is necessary but without it the plots don't work.
    par(new=TRUE)
    plot.new()

    ## Draw the transcript in the lower 1/3 of the plot area
    par(plt=c(plt[1], plt[2], plt[3], plt[3]+(plt[4]-plt[3])/3))
    par(usr=c(tmp[1], tmp[2], 0, 100))
    plot(transcript, xlim=xlim, units=units)

    ## restore the graphical environment
    par(usr=usr)
    par(plt=plt)
}


#gene="NR_029560"   # Mir150
#gene='NM_TEST'		# 'test'
gene='NM_016978'	# Oat '-'
gene='NM_144903'	# Aldob '-'
gene='NM_011434'	# 'Sod1'
gene='NM_011044'	# 'Pck1'
gene='NM_013541'	# 'Gstp1'
gene='NM_001005419'	# 'Ado'  single exon reverse
gene="NM_007409"   # ADH1

plot.new()

df = morris.getalignments("113010_A", gene)

##df <- data.frame(position=49557732, length=22)
##attr(df,"genome") <- "mm9"

kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
##kg <- data.frame(genome="mm9", name=c("gene1","gene2"), strand='+', txStart=1, txEnd=20,
##                 cdsStart=c(5,6), cdsEnd=c(8,9), exonCount=1, exonStarts="1,", exonEnds="20,",
##                 name2=c("test","test2"),stringsAsFactors = F)
rownames(kg) <- kg$name
##print(kg[gene,])

#df <- data.frame(position=c(49562310,49562310), length=22)
#attr(df,"genome") <- "mm9"

p = profile(df, kg[gene,])

plot(p, bias="middle", minlen=28, units="aa", xlim=c(15, 115))
plot(p, bias="middle", minlen=28, xlim=c(15, 115))


