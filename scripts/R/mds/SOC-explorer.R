
## read the SOC clinical outcomes spreadsheet.

## Identify interesting subsets of data to compare.

## create a command line for bam2depth tool using the names of bam files identified in the spreadhsheet

## call bam2depth to produce a tab-separated-values (.tsv) file suitable for use by the qarp functions.

library(qarp)

## read SOC clinical outcomes spreadsheet
prefix <- "/mnt/fred/"
if (!file.exists(paste0(prefix,"data"))) {
    prefix <- "~/"
}

df.survival <- read.csv(paste0(prefix, "data/SOC/2013-09-25_McIntoshSOC_Survival.csv"))
str(df.survival)
                      
levels(df.survival$Response)
levels(df.survival$Resection)

df.resistant <- subset(df.survival, Response == "Chemoresistant" & Resection == "Optimal")
df.sensative <- subset(df.survival, Response == "Chemosensitive" & Resection == "Optimal")

str(df.resistant)
str(df.sensative)

df.resistant <- df.resistant[1:10,]
df.sensative <- df.sensative[1:10,]

rfiles <-  paste0(gsub("_","-",df.resistant$Sample), ".bam")
sfiles <- paste0(gsub("_","-",df.sensative$Sample), ".bam")

# prepend the directry on each file
rfiles <- paste0(prefix, "data/SOC/", rfiles)
sfiles <- paste0(prefix, "data/SOC/", sfiles)

cmd <- 'python inst/python/bam2depth.py -s ~/src/samtools/samtools -g ACTB inst/python/hg19-refSeq-union.gtf'
cmd <- paste(cmd, paste(rfiles, collapse=" "), paste(sfiles, collapse=" "))
cmd <- paste(cmd, " > SOC-explorer.tsv")
print(cmd)
system(cmd)


