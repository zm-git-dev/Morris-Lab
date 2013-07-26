##Cairo(width = 640, height = 480, file="Rplot.pdf", type="pdf", pointsize=12,
##bg = "white", canvas = "white", units = "px", dpi = 72)
par(fig=c(0,0.50,0.1,0.9))

plot( x[,1], x[,2], cex=.5, log="xy",
panel.first = grid(4,4), cex.lab=.6,
xlab=paste0(descriptions[datasets[1],1], " (log)"),
ylab=paste0(descriptions[datasets[2],1], "(log)"))
par(fig=c(0.5,1,0.1,0.9), new=TRUE)

par(mgp=c(1,1,0))
plot( x[,1], x[,2], cex=.5, log="xy",
      panel.first = grid(4,4), cex.lab=.6, cex.axis=.5,
      xlab=paste0(descriptions[datasets[1],1], " (log)"),
      ylab=paste0(descriptions[datasets[2],1], "(log)"))
