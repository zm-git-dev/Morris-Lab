suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))

if (interactive()) {
} else {
    option_list <- list(
        make_option("--genome", default="unknown",
                    help="file containing name genes to monitor. [default \"%default\"]")
        )

    # get command line options, if help option encountered print help and exit,
    # otherwise if options not found on command line then set defaults,
    opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
}

if (length(opt$args) != 1) {
    stop("usage: union.R [options] <gene_model.gtf>")
}

gtfname <- opt$args[[1]]
genome <- opt$genome

gtf <- rtracklayer::import(gtfname, feature.type=c("exon"), asRangedData=FALSE, format="gtf")

genes <- unique(gtf$gene_id)
gl <- GRanges()
for (genename in genes) {
    print(genename)
    gene <- gtf[mcols(gtf)$gene_id==genename,]
    for (chrom in unique(seqnames(gene))) {
        cgene <- gene[seqnames(gene)==chrom,]
        for (s in unique(strand(cgene))) {
            sgene <- cgene[strand(cgene)==s,]
            rgene <- reduce(sgene)
            mcols(rgene) <- list(gene_id=genename,
                         transcript_id=paste0(unique(sgene$transcript_id), collapse=":"),
                         type="exon", exon_number=seq(1:length(rgene)))
            gl <- append(gl, rgene)
        }
    }
}

export(gl,  paste0(genome, "_union.gtf"), format="gtf", source=paste0(genome, "_union"))

genename <- "LACTB"
