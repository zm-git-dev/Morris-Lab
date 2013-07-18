
##library(grid)


## define S3 class 'transcript' to hold instances gene transcript info.
##
## generally not used directly but it used as part of the 'profile' class.
##
## Typical vignette:
##    gene ="NM_009654"   ## Albumin
##    df = morris.getalignments("113010_A", gene)
##    kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
##    p = profile(df, kg[gene,])
##    print(p)
##    plot(p, minlen=28)
##
## This awesome flavor of OOP in R using closures is taken from a
## stackoverflow posting I stumbled across:
## http://stackoverflow.com/a/15245568/1135316
## I've since realize that this is just the S3 style of classes in R.
## Add a few lines and it is now an S3 style class.
transcript = function(gdata) {
    starts <- as.numeric(strsplit(gdata$exonStarts,",")[[1]])
    ends <- as.numeric(strsplit(gdata$exonEnds,",")[[1]])
    elens <- ends - starts
    transcriptLength <- sum(elens)

    cdsStart = function() rpos(if (gdata$strand == '+') { gdata$cdsStart } else { gdata$cdsEnd })
    cdsEnd = function() rpos(if (gdata$strand == '+') { gdata$cdsEnd-1 } else { gdata$cdsStart+1 })
    cdsLength = function() cdsEnd()-cdsStart()+1
    peptideLength = function() cdsLength()/3
    txStart = function()  rpos(if (gdata$strand == '+') { gdata$txStart } else { gdata$txEnd })
    txEnd = function() rpos(if (gdata$strand == '+') { gdata$txEnd-1 } else { gdata$txStart+1 })
    txLength = function() txEnd()-txStart()+1
    isCoding = function() (gdata$cdsEnd != gdata$cdsStart)
    name = function() (gdata$name)
    name2 = function() (gdata$name2)
    cpos = function(offset) {
        ## Give a an offset from the beginning of the transcript,
        ## return an absolute chromosome position.
        pos = offset
        pos = pos - 1
        if (gdata$strand == '-') {
            pos = transcriptLength + pos - 1
        }
        for (i in seq(length(elens))) {
            elen = elens[i]
            if (pos - elen < 0)
              break
            pos = pos - elen
        }
        pos = pos + starts[i]
        return(pos)
    }
    
    rpos = function(pos) {
        ## Given a position on a chromosome, calculate the position
        ## relative to the start of transcript.
        ##
        ## relative position is equal to the sum of all exons that end
        ## before the position, PLUS the beginning portion of the exon
        ## that contains the position.

        offset = 0
        ## there are two primary situations:
        ## 1) the position is to the right of the start of the transcript
        ## 2) the position is to the left of the start of the transcript
        if (pos >= starts[1]) {
            ## positioned to the right of the beginning of the gene
            if (any(ends<=pos)) {
                offset = sum( elens[ends<=pos] )
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
    
    stopifnot(transcriptLength == txLength())
    stopifnot(txStart() == 1)
    
    exported = list(
      txLength=txLength,
      cdsLength=cdsLength,
      cdsStart=cdsStart,
      cdsEnd=cdsEnd,
      txStart=txStart,
      txEnd=txEnd,
      isCoding=isCoding,
      name=name,
      name2=name2,
      rpos=rpos,
      cpos=cpos,
      peptideLength=peptideLength
      )
    class(exported) <- "transcript"
    invisible(exported)
}

## specialization of generic print function for instances of 'transcript' class.
print.transcript <- function(this) {
    cat(sprintf("gene %s\n", this$name2()))
    cat(sprintf("   tx %d-%d  (%d)\n", this$txStart(), this$txEnd(), this$txLength()))
    cat(sprintf("   cds %d-%d  (%d)\n", this$cdsStart(), this$cdsEnd(), this$cdsLength()))
}

## specialization of generic plot function for instances of 'transcript' class.
## produces a plot that looks something like this, 
##
##        -------|||||||||||||--------
##
## where the thin part (dashes) represents the 5'- and 3'-UTR, and the
## thick part (vertical bars) represents the coding region of the
## transcript.
##
plot.transcript <- function(self, xlim=NULL, units="nucleotide") {
    if (units == "aa") {
        if (is.null(xlim)) 
          xlim <- c(self$cdsStart(), self$cdsEnd())
    } else {
        if (is.null(xlim))
          xlim <- c(1, self$txLength())
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
        usr = par()$usr
        clip(xlim[1], xlim[2], 0, 100)
        segments(xlim[1], 100-cdsHeight,  xlim[2],  100-cdsHeight, lwd=5, lend="butt")
        rect(cdslim[1], 100-2*cdsHeight, cdslim[2], 100, col='blue')
        do.call("clip", as.list(usr))  # reset to plot region
    } else {
        ## There is no coding region - this is a non-coding gene.
        ## Label the endpoints of the transcript and label the gene in the middle.
        text((xlim[1]+xlim[2])/2, 100-2*cdsHeight-5,
             adj=c(.5,1.0), self$name2())
        text(xlim[1], 100-2*cdsHeight-5, adj=c(.5,1.0), as.character(xlim[1]))
        text(xlim[2], 100-2*cdsHeight-5, adj=c(.5,1.0), as.character(xlim[2]))
        usr = par()$usr
        clip(xlim[1], xlim[2], 0, 100)
        segments(xlim[1], 100-cdsHeight,  xlim[2],  100-cdsHeight, lwd=5, lend="butt")
        do.call("clip", as.list(usr))  # reset to plot region

    }
}

## define S3 class 'profile' to hold instances of ribosome profile data.
## A 'profile' is the set of ribosome positions for a gene in an experiment.
## One profile instance would hold data for one gene in one experiment.
##
## Typical vignette:
##    gene ="NM_009654"   ## Albumin
##    df = morris.getalignments("113010_A", gene)
##    kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
##    p = profile(df, kg[gene,])
##    print(p)
##    plot(p, minlen=28)
##
profile <- function(df, gene.data) {
    xsript <- transcript(gene.data)
    transcript <- function() xsript
    alignments <- function() df
    rpositions <- function() {   # relative positions w.r.t. start of transcript
        rpositions = sapply(X=df$position, FUN=xsript$rpos)
    }
    plotpositions <- function(xlim=NULL, units="nucleotide", bias="middle", minlen=0) {
        ## filtered plot positions
        df <- alignments()
        transcript <- transcript()
        
        ## Calculate the relative position of each alignment w.r.t. start of transcript.
        ## position may be calculated with respect to the middle of the read or
        ## the 5'-end
        if (bias == "middle") 
          df$rpositions = df$rpositions + df$length/2.0

        ## if the user is zooming in on a portion of the transcript, set the
        ## limits of the horizontal axis accordingly.   Otherwise the limits
        ## are determined by the length of the transcript.
        if (units == "aa") {
            ## convert nucleotide positions to codon position
            df$rposition = round((df$rposition - transcript$cdsStart())/3)
        } 
    
        ## discard any alignments that fall outside the horizontal limits.
        ## ( why?  won't the plot just truncate any data outside the range? )
        if (!is.null(xlim)) {
            is.between <- function(x,lim) {
                (x > lim[1]) & (x < lim[2])
            }
            df <- df[is.between(df$rposition,xlim),]
        }

        ## if the user specified a minimum read length, discard shorter reads.
        ##print (df[df$len < minlen,])
        if (minlen < 0) 
          df = df[df$len < -minlen,]
        else
          df = df[df$len >= minlen,]
        
        return(df)
    }

    ## precalculate the relative ribosome positions on the transcript.
    df$rpositions = rpositions()
    
    ## make a list of methods exported from instances of this class.
    exported = list(
      transcript=transcript,
      alignments=alignments,
      rpositions=rpositions,
      plotpositions=plotpositions )

    ## turn into an S3 class by setting the class attribute.
    class(exported) <- "profile"

    ## return a list of closures, with some local variables already bound.
    ## that's all an S3 class is - a list of closures.
    invisible(exported)
}

# Summary method
# Should return a summary object (and not print it!)
# Need a unique class for it ("fsummary")
summary.profile <- function(self,...) {
    res <- summary(self$alignments())
    res
}


## specialization of generic print function for instances of 'profile' class.
print.profile <- function(this) {
    print(this$transcript())
}

## specialization of generic plot function for instances of 'profile' class.
plot.profile <- function(self, xlim=NULL, units="nucleotide", bias="middle", minlen=0, identify=FALSE, ...) {
    usr = par()$usr
    plt = par()$plt

    transcript <- self$transcript()
    ## if the user is zooming in on a portion of the transcript, set the
    ## limits of the horizontal axis accordingly.   Otherwise the limits
    ## are determined by the length of the transcript.
    ## Calculate the histogram bins and maximum count across all bins.
    ##
    complete = FALSE
    if (units == "aa") {
        if (is.null(xlim)) {
            complete = TRUE
            xlim = c(1, transcript$peptideLength())
        } else if (length(xlim) == 1) {
            xlim = c(xlim[1], transcript$peptideLength())
        }
        breaks=seq(xlim[1], xlim[2], 1)
    } else {
        if (is.null(xlim)) {
            complete = TRUE
            xlim <- c(1, transcript$txLength())
        } else if (length(xlim) == 1) {
            xlim = c(xlim[1], transcript$txLength())
        }
        ## for small transcripts with few reads, use a bin size of one.
        ## otherwise use a bin size of three
        if ((xlim[2] - xlim[1]) > 100) 
          breaks=seq(xlim[1],xlim[2]+1, 3)
        else
          breaks=seq(xlim[1],xlim[2], 1)
    }

    ## calculate the ribosome positions relative to start of transcript.
    df <- self$plotpositions(xlim=xlim, units=units, bias=bias, minlen=minlen)
    
    ## Draw the histogram in the top 2/3 of the plot area.
    par(plt=c(plt[1], plt[2], plt[3]+(plt[4]-plt[3])/3, plt[4]))

    histdata <- hist(df$rposition, breaks=breaks, plot=F)
    maxy <- max(histdata$count)

    plot(histdata$mids[histdata$counts != 0],histdata$counts[histdata$counts != 0], type='h',
         xlim=xlim, ## ylim=c(0,1.25*maxy), 
         lwd=3, lend=2,ylab='',xlab='',main='',xaxt='n',bty='n', ...)

    # The usr coordinates appear to change after plot(), so fetch the coordinates
    # again before positioning the labels.
    tmp = par("usr")
    text(xlim[1], tmp[4], attr(df, "dataset"), adj=c(0,1.2), new=TRUE)
    text(xlim[2], tmp[4], paste0(nrow(df), " total reads"), adj = c( 1, 1.2 ), new=TRUE)

    ## start a new plot without clearing the current one.
    ## I don't understand why this is necessary but without it the plots don't work.
    par(new=TRUE)
    plot.new()

    ## Draw the transcript in the lower 1/3 of the plot area
    ## plt(x1, x2, y1, y2)
    par(plt=c(plt[1], plt[2], plt[3], plt[3]+(plt[4]-plt[3])/3))

    ## usr(x1, x2, y1, y2)
    par(usr=c(tmp[1], tmp[2], 0, 100))
    plot(transcript, xlim=xlim, units=units)

    legend = ""
    legend=paste(legend,"units=",paste0("\"", units, "\""))
    if (complete) {
        legend=paste(legend,"whole transcript")
    } else {
        legend=paste(legend,",xlim=",xlim[1], ",", xlim[2])
    }
    legend=paste(legend,"bias=",paste0("\"", bias, "\""))
    legend=paste(legend,"minlen=",minlen)
    
    text(tmp[1], 0, legend, adj=c(0,-1), col="grey56")
    
    par(plt=c(plt[1], plt[2], plt[3]+(plt[4]-plt[3])/3, plt[4]))
    par(usr=tmp)

    if (identify == TRUE) {
        ## A function to use identify to select points, and overplot the
        ## points with another symbol as they are selected
        identifyPch <- function(x, y=NULL, n=length(x), pch=19, ...)
          {
              xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
              sel <- rep(FALSE, length(x)); res <- integer(0)
          }

        identifyPch(histdata$mids[histdata$counts != 0],histdata$counts[histdata$counts != 0])
        
        xy <- xy.coords(histdata$mids[histdata$counts != 0],histdata$counts[histdata$counts != 0])
        x <- xy$x; y <- xy$y
        
        repeat {
            ans <- identify(x, y, n=1, plot=FALSE, ...)
            if(!length(ans)) break
            print(paste0("(",x[ans], ",", y[ans],")  ", self$transcript()$cpos(x[ans])))
        }
    }
    
    ## repeat {
    ##     coord = locator(1)
    ##     if(is.null(coord)) 
    ##         break
    ##     print(paste0(round(coord$x), ",", round(coord$y)))
    ## }

    ## restore the graphical environment
    par(usr=usr)
    par(plt=plt)
    invisible(df)
}



