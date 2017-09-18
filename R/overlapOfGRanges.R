overlapOfGRanges <- function(GRlist, plot=TRUE) {
    if(!is(GRlist, 'list')) stop('GRlist has to be a list ...')
    if(!all(sapply(GRlist, function(x) is(x, 'GRanges'))))
        stop('GRlist has to be a list of GRanges ...')
    if(!is.logical(plot))
        stop('plot has to be of class logical ...')

    grL <- length(GRlist)
    overlapMatrix <- matrix(100, grL, grL)
    rownames(overlapMatrix) <- names(GRlist)
    colnames(overlapMatrix) <- names(GRlist)
    for(i in 1:grL) {
        for(j in (1:grL)[-i]) {
            res <- countOverlaps(query= GRlist[[i]], subject= GRlist[[j]])
            count <- length(which(res > 0))
            overlapMatrix[i,j] <- 100 * count / length(GRlist[[i]])
        }
    }
    if(plot) {
        cpWB <- colorRampPalette(c('white','beige'))
        cpBR <- colorRampPalette(c('beige','red'))
        cols <- c(cpWB(50), cpBR(51)[-1])
        heatmap.2(overlapMatrix, cexRow=.7, cexCol=.7, col=cols, Rowv=NULL,
                  Colv=NULL, cellnote=signif(overlapMatrix, 2),
                  dendrogram='none', trace='none', density.info='none', notecol='black')
    }

    invisible(overlapMatrix)
}
