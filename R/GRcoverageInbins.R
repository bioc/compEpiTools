setGeneric('GRcoverageInbins', function(Object, bam, Nnorm=TRUE,
                                        Snorm=FALSE, Nbins) standardGeneric('GRcoverageInbins'))
setMethod('GRcoverageInbins','GRanges', function(Object, bam,
                                                 Nnorm=TRUE, Snorm=FALSE, Nbins) {
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
    widths <- width(Object)
    startPosMat= matrix(NA, length(Object), Nbins)
    for(i in 1:length(Object)) startPosMat[i,]= round(seq(1,widths[i], length.out=Nbins+1))[-(Nbins+1)]
    endPosMat= cbind(startPosMat[,-1]-1, as.numeric(widths))
    binsizeMat= endPosMat-startPosMat+1
    for(bin in 1:Nbins) {
        for(gr in 1:length(Object)) coverageMat[gr,bin] <-
            sum(covList[[gr]][startPosMat[gr,bin]:endPosMat[gr,bin]])
        if(Snorm) {
            coverageMat[,bin] <- coverageMat[,bin] / binsizeMat[,bin]
        }
    }
    coverageMat[widths<Nbins,]<- NA
    return(coverageMat)
})

