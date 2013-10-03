
##  Time-stamp: <2013-08-26 14:56:52 chris>

## Create a scatter plot of gene expression for prostate experiments.
##   -- eliminate all Riken genes
##   -- eliminate non-coding genes.
##   -- eliminate lowest expression TR- control sets. (<500K total reads)
##   -- eliminate the 3 least developed tumors from the TR+ sample set
##   -- highlights the fibroblast, epithelial, and tumor-associated genes.
##
## This script will draw a scatter plot of gene expression levels and
## will overlay a regession line that passes through the sibroblasts
## genes.  Then, for 10 regions along the x-axis, it draws two
## flanking lines that indicate the mean of the absolute difference
## from each point in the region to the regression line.
##

 

source("../morrislib.R")
library(ggplot2)
library(Biobase)
library(Hmisc)
require(genefilter)

## These are color-blind-friendly palettes, one with gray, and one with black.
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tumorSpecific <- c( "Tubb5", "Hist2h2bb", "Nov", "Chgb", "Hist1h1a" )
epithelial <- c("Tgm4", "Pbsn", "Sbp", "Rnase1")
fibroblast <- c("Vim", "Itgb1", "Itga1", "Col1a1", "Col1a2")

# divide the data into two sets - control and treated...
##treated <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+/TR+")
treated <- c("103112_A", "030713_B", "041713_B")
control <- morris.datasets(organism="Mouse", tissue="Prostate", genotype="Ribo+/Col+")
datasets <- c(treated,control)
df.raw = morris.genecounts(datasets, group=NULL)
genome = attr(df.raw, "genome")

kg <- morris.getknowngenes(genome)

## eliminate lowest expression TR- control sets. (<600K total reads)
control <- control[colSums(df.raw[control], na.rm = TRUE) > 600000]

df <- df.raw[, c(control, treated)]

## Add a pseudocount to all missing data so logarithms taken later
## will be well-defined.
df[is.na(df)] <- 1

## toss out non-coding genes
df <- df[grep("NR_", rownames(df), invert=TRUE),]

## toss RIKEN genes.
df <- df[grep("Rik", kg$name2[match(rownames(df), kg$name)], invert=TRUE), ]

## toss any gene that does not have a corresponding common name
df <- df[!is.na(kg$name2[match(rownames(df), kg$name)]), ]

## scale all gene counts to the count of Col1a2
Col1a2 = with(kg, name[match("Col1a2", name2)])
n = as.numeric(df[Col1a2, control])
df[control] = sweep(df[control], 2, n, '/')
df[control] = sweep(df[control], 2, mean(n), '*')
n = as.numeric(df[Col1a2, treated])
df[treated] = sweep(df[treated], 2, n, '/')
df[treated] = sweep(df[treated], 2, mean(n), '*')

##df <- df[apply(df, 1, function(row) (!all(row < 40))),]
## eliminate genes for which the average control and average treated are BOTH below 40 reads
df <- df[!(rowMeans(df[control]) < 40 & rowMeans(df[treated]) < 40),]

## Normalize to mapped reads per million mapped reads
attr(df, "genome") = genome
df.norm = morris.normalize(df[treated], normalization="rpm")
df.norm[control] = morris.normalize(df[control], normalization="rpm")

names(epithelial) = with(kg, name[match(epithelial, name2)])
names(fibroblast) = with(kg, name[match(fibroblast, name2)])
names(tumorSpecific) = with(kg, name[match(tumorSpecific, name2)])

df.norm$tissue <- "0"
df.norm[ intersect(names(epithelial), rownames(df.norm)), "tissue" ] <- "epithelial"
df.norm[ intersect(names(fibroblast), rownames(df.norm)), "tissue" ] <- "fibroblast"
df.norm[ intersect(names(tumorSpecific), rownames(df.norm)), "tissue" ] <- "tumor"
df.norm$tissue = as.factor(df.norm$tissue)

df.norm$m1 = log2(rowMeans(df.norm[,control]))
df.norm$m2 = log2(rowMeans(df.norm[,treated]))

shapes <- c("0"=20, "epithelial" = 20, "fibroblast" = 15, "tumor"=17)
sizes <- c("0"=2, "epithelial" = 4, "fibroblast" = 4, "tumor"=4)
cbbPalette <- c("lightgrey", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## calculate a regression line through the fibroblast genes.
model = lm( m2 ~ m1, data=df.norm[names(fibroblast),])

## calculate error estimates for every point in m1
fun = function(x,y,model) {
  abs(y - predict(model, data.frame(m1=x)))
}
eest = mapply(fun, df.norm$m1, df.norm$m2, MoreArgs = list(model=model), SIMPLIFY=TRUE)

## Divide m1 into bins
bpts = pretty(df.norm$m1, n=10)
INDEX = cut(df.norm$m1, bpts, include.lowest = T)

## caclulate the average x value for each of the bins
x = tapply(df.norm$m1, INDEX, mean)

## calculate the average error for each of the bins
err = tapply(eest, INDEX, mean)

# calcuate predicted points along the line for each of the x values
y = sapply(x, function(x) predict(model, data.frame(m1=x)))

gg <- ggplot(df.norm[with(df.norm, order(tissue)), ], aes(name="", x=m1, y=m2))
gg <- gg + theme_bw()
gg <- gg + theme(legend.position="top",legend.title=element_blank(), legend.text = element_text(colour="blue", size = 14, face = "bold"))
gg <- gg + ylab("Tramp+ (log2 rpm)")
gg <- gg + xlab("Control (log2 rpm)")
gg <- gg + geom_point(aes(shape = tissue, color = tissue, size=tissue)) 
gg <- gg + scale_size_manual(name="", values=sizes, breaks=c("epithelial", "fibroblast",  "tumor")) 
gg <- gg + scale_color_manual(name="", values=cbbPalette, breaks=c("epithelial", "fibroblast",  "tumor"))
gg <- gg + scale_shape_manual(name="", values=shapes, breaks=c("epithelial", "fibroblast",  "tumor"))
gg <- gg + geom_point(data=df.norm[names(fibroblast),], size=4, col="red")
gg <- gg + geom_smooth(method=lm, data=df.norm[names(fibroblast),], se=FALSE, fullrange=TRUE)
gg <- gg + geom_line(data=data.frame(x=x, y=y+err), aes(x=x, y=y))
gg <- gg + geom_line(data=data.frame(x=x, y=y-err), aes(x=x, y=y))
print(gg)

## Sys.sleep(0) 

## repeat {
##     ans <- identify(m1, m2, n=1, plot=FALSE)
##     if(!length(ans)) break
##     print(paste0("(",m1[ans], ",", m2[ans],")  ",rownames(df.norm)[ans]," ", with(kg, name2[match(rownames(df.norm)[ans], name)])))
## }

## Sys.sleep(0) 
## s1 = rowSds(df.raw[,control])
## s2 = rowSds(df.raw[,treated])
## ttest=(df.norm$m2-df.norm$m1)/sqrt(s1^2/length(control)+s2^2/length(treated))
## hist(ttest,nclass=100)

## Sys.sleep(0) 
## ## Calculate the two-tailed T-test statistic
## ## probabiliy of seeing the a value greater than the measured value on 
## ## both the positive and negative extremese.
## pval=2*(1-pt(abs(ttest),4))
## hist(pval, nclass=100)
## Sys.sleep(0) 


