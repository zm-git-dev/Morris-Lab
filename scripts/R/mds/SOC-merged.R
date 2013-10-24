library(qarp)

data <- qarp.read.cachedtsv("/Volumes/homes/data/SOC-merged.tsv")
control <- c(
    "Benign_1312_501369.depth", "Benign_1675_224804.depth", "Benign_441_206864.depth",
    "Benign_2638_223151.depth", "Benign_7012_315111.depth", "Benign_7609_361548.depth",
    "Benign_4764_251506.depth", "Benign_1069_501203.depth", "Benign_116_100831.depth",
    "Benign_2186_509428.depth")

treated <- c(
    "SOC_7637_361542.depth", "SOC_7777_371281.depth", "SOC_12523_494920.depth",
    "SOC_849_206653.depth", "SOC_9547_467919.depth", "SOC_6208_299803.depth",
    "SOC_5991_294171.depth", "SOC_2745_226696.depth", "SOC_5959_278682.depth",
    "SOC_13451_492771.depth")

datasets = c(control = control, treated = treated)

# eliminate extra columns that we don't need.
data = data[, c("symbol", "transPos", "exonNumber", datasets)]
data$symbol <- as.factor(data$symbol)

aVec <- as.factor(sapply(datasets %in% control, function(x) if (x) {"control"} else {"treated"}))
names(aVec) <- datasets


pvalues = qarp(data, aVec)


mat <- qarp.matrix(data, "LCN2", aVec, meancenter=FALSE)
qarp.plotProfile(mat, aVec)
#qarp.plotMDS(mat, aVec)

mat <- qarp.matrix(data, "ACTB", aVec, meancenter=TRUE)
#qarp.plotMDS(mat, aVec)

pvalues <- qarp(data, aVec, genelist=c("ACTB"))
print(pvalues)


df = qarp(data, aVec)
print(str(df))
#hist(df$pvalues, n=50, main=attr(data, "filename"))
df <- df[order(df$pvalues),,drop=FALSE]
print(df)

