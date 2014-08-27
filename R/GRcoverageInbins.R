setGeneric('GRcoverageInbins', function(Object, bam, Nnorm=TRUE,
                                        Snorm=TRUE, Nbins) standardGeneric('GRcoverageInbins'))
setMethod('GRcoverageInbins','GRanges', function(Object, bam,
                                                 Nnorm=TRUE, Snorm=TRUE, Nbins) {
    if(!is.character(bam))
        stop('bam has to be a file path of class character ..')
    if(!is.logical(Nnorm) || !is.logical(Snorm))
        stop('Nnorm and Snorm have to be of class logical ..')
    if(!is.numeric(Nbins))
        stop('Nbins has to be of class numeric ..')
    if(Nbins>min(width(Object)))
        warning('Nbins is greater than the width of some ranges, this will result in NAs for those ranges ...')

    covList <- GRbaseCoverage(Object=Object, bam=bam, Nnorm=Nnorm)

    coverageMat <- matrix(NA, length(Object), Nbins)
    binsize <- round(width(Object) / Nbins)
    widths <- width(Object)
    for(bin in 1:Nbins) {
        startPos <- binsize * (bin - 1) + 1
        if(bin == Nbins) endPos <- widths
        else endPos <- startPos + binsize - 1
        for(gr in 1:length(Object)) coverageMat[gr,bin] <-
            sum(covList[[gr]][startPos[gr]:endPos[gr]])
        if(Snorm) {
            if(bin == Nbins) coverageMat[,bin] <-
                coverageMat[,bin] / (widths- binsize*(Nbins-1))
            else coverageMat[,bin] <- coverageMat[,bin] / binsize
        }
    }
    coverageMat[Nbins>widths]<- NA
    return(coverageMat)
})

