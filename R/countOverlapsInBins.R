setGeneric('countOverlapsInBins', function(query, subject, nbins)
           standardGeneric('countOverlapsInBins'))
setMethod('countOverlapsInBins','GRanges',
          function(query, subject, nbins) {
    if(!is(subject,'GRanges'))
        stop('subject has to be of class GRanges ..')
    if(!is.numeric(nbins))
        stop('nbins has to be of class numeric ..')
    countMat<- matrix(0, length(query), nbins)
    binsize<- floor(width(query)/nbins)
    starts<- start(query)
    for(bin in 1:nbins) {
        startPos<- starts+ binsize*(bin-1)
        if(bin == nbins) endPos<- end(query)
        else endPos<- startPos + binsize - 1
        queryBin<- GRanges(
            seqnames= seqnames(query),
            ranges=IRanges(start=startPos, end=endPos))
        countMat[, bin]<- countOverlaps(queryBin, subject)
    }
    countMat[countMat > 1]<- 1

    return(countMat)
})

