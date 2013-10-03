
## plot the number of reads discovered for each experiment using three
## different criteria.

df = read.table("table.txt", row.names=1, stringsAsFactors = FALSE)
                
                
# Grouped Bar Plot

barplot(t(as.matrix(df)), main="Read count by experiment and discriminant",
        xlab="Experiment", 
        legend = rownames(counts), beside=TRUE, col=heat.colors(3))

# Place the legend at the top-left corner with no frame  
# using rainbow colors
legend("topleft", c("no mismatch","1 mismatch","first seven"), cex=0.6, 
       bty="n", fill=heat.colors(3));