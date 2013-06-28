

group=NULL
treated <- c("061113_B")
control <- c("061113_D")
datasets <- c(control)

datasets <- c('061113_A', '061113_B','061113_C','061113_D')

foo = function() {
    ds <- morris.genecounts(datasets, group=group)
    genome <- attr(ds, "genome")
    stopifnot(!is.null(genome))
}
debug("foo")
foo()
rpkm <- morris.normalize(df, normalization="rpkm", group=group)
rownames(rpkm)[!complete.cases(rpkm)]
