

samples <- read.csv("SOC_pvalues.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, comment.char="#")
# some genes will have '*' in the name if some experiments have no reads for that gene
# strip and ignore this flag for now.
samples[,1] <- sub("*", "", samples[,1], fixed=TRUE)
samples[,2] = as.numeric(samples[,2])
samples[,3] = as.numeric(samples[,3])
samples[,4] = as.numeric(samples[,4])

png("pvalue.qqplots.png", width=680, height=680)
panel.qqplot <- function(x, y, ...) {
  lines(qqplot(x,y,plot.it=FALSE), lwd=3, ...)
}
pairs(samples[,2:4], panel=panel.qqplot,, main="qqplots of pvalue distributions\nfor various normalization methods")
dev.off()


png("pvalue.histplots.png", width=680, height=680)
par(mfrow=c(2,2), oma=c(2.5, 3, 2.5, 2.5), mar=c( .1, .1, 2.5, .1), cex=1, las=1)
hist(samples[,2], yaxt='n', xlab="pvalue", main="")
mtext(sprintf("Not normalized"), side=3, line=-2, cex=1.3)

hist(samples[,3], yaxt='n', xlab="pvalue", main="")
mtext(sprintf("RPM normalized"), side=3, line=-2, cex=1.3)

hist(samples[,4], yaxt='n', xlab="pvalue", main="")
mtext(sprintf("Area normalized"), side=3, line=-2, cex=1.3)

mtext(sprintf("Distribution of pvalues\n%d genes",  nrow(samples)), side=3, line=-1, cex=1.3, outer=TRUE)
dev.off()
