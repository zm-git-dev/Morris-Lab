library(qarp)

# different prefixes allow me to operate on both the clister servers and on my laptop.
prefix <- "/mnt/fred/"
if (!file.exists(paste0(prefix,"data"))) {
    prefix <- "~/"
}

data <- qarp.read.cachedtsv("SOC-explorer.tsv")

## read SOC clinical outcomes spreadsheet
df.survival <- read.csv(paste0(prefix, "data/SOC/2013-09-25_McIntoshSOC_Survival.csv"))

## Gather patient records that differ by chemoresistance
df.resistant <- subset(df.survival, Response == "Chemoresistant" & Resection == "Optimal")
df.sensative <- subset(df.survival, Response == "Chemosensitive" & Resection == "Optimal")

## only use the first 10 record from each set (to reduce time in analysis)
df.resistant <- df.resistant[1:10,]
df.sensative <- df.sensative[1:10,]


datasets = c(control = df.resistant, treated = df.sensative)

total.aligned <- c("SOC.2745.226696"=14612031, "SOC.849.206653"=14494734)


# eliminate extra columns that we don't need.
data = data[, c("symbol", "transPos", "exonNumber", datasets)]
data$symbol <- as.factor(data$symbol)

aVec <- as.factor(sapply(datasets %in% control, function(x) if (x) {"control"} else {"treated"}))
names(aVec) <- datasets


pvalues = qarp(data, aVec)

mat <- qarp.matrix(data, "ACTB", aVec)
print(apply(mat,1,sum)/1e6)

# scale by the total number of reads in the gene.
# mat <- t(scale(t(mat),center=FALSE,scale=apply(mat,1,sum)/1e6))

# scale by the total number of aligned reads in the entire experiment.
totals <- total.aligned[rownames(mat)]/1e6
mat <- t(scale(t(mat),center=FALSE,scale=totals))

print(str(mat))
qarp.plotProfile(mat, aVec, invert=FALSE)
#qarp.plotMDS(mat, aVec)

pvalues <- qarp(data, aVec, genelist=c("ACTB"))
print(pvalues)


df = qarp(data, aVec)
print(str(df))
#hist(df$pvalues, n=50, main=attr(data, "filename"))
df <- df[order(df$pvalues),,drop=FALSE]
print(df)

