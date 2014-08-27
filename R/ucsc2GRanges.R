ucsc2GRanges <- function(ucscPositions) {
    if(!is.character(ucscPositions))
        stop('ucscPositions has to be of class character')

    xsplit <- strsplit(ucscPositions, split=':')
    xchr <- sapply(xsplit, function(x) x[1])
    xse <- sapply(xsplit, function(x) x[2])
    xse <- strsplit(xse, split='-')
    xs <- as.numeric(sapply(xse, function(x) x[1]))
    xe <- as.numeric(sapply(xse, function(x) x[2]))
    gr <- GRanges(seqnames=xchr, IRanges(start=xs, end=xe))

    return(gr)
}
