
library(grid)

.pardefault <- par(no.readonly = T)


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
        query <- paste0("select position,length from ",
                        "new_alignments_tbl a join datasets_tbl d on a.dataset_id=d.id ",
                        "where d.name like '", dataset, "%' and a.feature like '", gene, "'")
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
    return(exported)
}

print.transcript <- function(object) {
    cat("transcript ", object$name2(), ", ", object$length(), "\n")
}

plot.transcript <- function(obj) {
    ## draw the transcript in the third frame down.
    plot(c(0,  obj$length()), c(0,100), type="n", axes=F, xlab='', ylab='', bty="n", new=T)

    cdsHeight = 5   # height of coding region, expressed as percentage of plot height
    segments(0, 100-cdsHeight,  obj$length(),  100-cdsHeight, lwd=5, lend=2)
    if (obj$isCoding()) {
        ## Label the endpoints of the coding region and place the gene
        ## name in the middle of the coding region.
    
        rect(obj$cdsStart(), 100-2*cdsHeight, obj$cdsEnd(), 100, col='blue')
        text((obj$cdsStart() + obj$cdsEnd())/2, 100-2*cdsHeight-5,
             adj=c(.5,0.5), obj$name2())
        text(obj$cdsStart(), 100-2*cdsHeight-5, adj=c(.5,0.5), as.character(obj$cdsStart()))
        text(obj$cdsEnd(), 100-2*cdsHeight-5, adj=c(.5,0.5), as.character(obj$cdsEnd()))
    } else {
        ## There is no coding region - this is a non-coding gene.
        ## Label the endpoints of the transcript and label the gene in the middle.
        text( obj$length()/2, 100-2*cdsHeight-5,
             adj=c(.5,0.5), obj$name2())
        text(0, 100-2*cdsHeight-5, adj=c(.5,0.5), as.character(0))
        text( obj$length(), 100-2*cdsHeight-5, adj=c(.5,0.5), as.character( obj$length() ))
    }
}


#gene="NR_029560"   # Mir150
gene="NM_007409"   # ADH1
gene='NM_001005419'	# 'Ado'  single exon reverse
gene='NM_013541'	# 'Gstp1'
gene='NM_016978'	# 'Oat'
gene='NM_011434'	# 'Sod1'
gene='NM_144903'	# Aldob
gene='NM_011044'	# 'Pck1'

# read all the alignments from 'dataset' that align to 'gene'
# so these will be multiple alignments, but only those that align to
# a single gene.
#
df = morris.getalignments("113010_A", gene)

# read all known genes
kg <- morris.getknowngenes(attr(df, "genome"), gene=gene, group=NULL)
rownames(kg) <- kg$name

gobj = transcript(kg[gene,])

## For each alignment from 'gene', add a column for position within that transcript.
df$transcriptPosition = sapply(X=df$position, FUN=gobj$rpos)

par(.pardefault)

#par(fig=c(0,0.2,0,0.2), new=FALSE)
plotIDs <- matrix(c(1:4), 4, 1, byrow=T)
layout(plotIDs, widths = c(1), heights = c(0.5,1,1,0.5))
par(mai=c(0, 0.5, 0, 0))
frame()  ## skip the top-most frame
## draw a histogram in the second frame down.

## Calculate the histogram bins and maximum count across all bins.
##
histdata = hist(df$transcriptPosition,breaks= gobj$length()/3,plot=F)
maxy=max(histdata$count[histdata$count != 0])
par(usr = c(0,  gobj$length(), 0, maxy) )

plot(histdata$mid[histdata$count != 0],histdata$count[histdata$count != 0], type='h',
      xlim=c(0, gobj$length()), ylim=c(0,maxy+10), 
      lwd=3, lend=2,ylab='',xlab='',main='',xaxt='n',bty='n')


# Put the total number of reads in the upper-right corner of the histogram.
usr <- par( "usr" )
text( usr[ 2 ], usr[ 4 ], paste0(nrow(df), " total reads"),    adj = c( 1.2, 5 ))
 
plot(gobj)

