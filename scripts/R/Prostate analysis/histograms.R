resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}


datasets <- c("103112_B_MM1",
              "030713_A_MM1",
              "032513_B_MM1",
              "103112_A_MM1",
              "030713_B_MM1",
              "032513_A_MM1")
df = morris.genecounts(datasets, group=NULL)
par(mfrow=c(2,3))
for (n in names(e)) {
  hist(e[,n], main=n, xlab="coverage depth")
  
}
par(mfrow=c(1,1))

