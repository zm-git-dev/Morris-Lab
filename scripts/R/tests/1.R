
## test script for morris.getGenome.....

test1 <- function() {
    datasets=c("113010_A", "073012_B")
    gs <- morris.getGenome(datasets)
    print(gs)
    return (length(unique(gs)) != 1)
}
test2 <- function() {
    datasets=c("NONE")
    gs <- morris.getGenome(datasets)
    print(gs)
    return (is.null(gs))
}
test3 <- function() {
    datasets=c("113010_A")
    gs <- morris.getGenome(datasets)
    print(gs)
    return (length(unique(gs)) == 1)
}
test4 <- function() {
    datasets=c("NONE", "113010_A")
    gs <- morris.getGenome(datasets)
    print(gs)
    return (length(unique(gs)) == 1)
}

test.genome <- function()  {
    checkTrue(test1())
    checkTrue(test2())
    checkTrue(test3())
    checkTrue(test4())
}
