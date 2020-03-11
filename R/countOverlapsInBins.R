setGeneric('countOverlapsInBins', function(query, subject, nbins,strandspecific=FALSE)
           standardGeneric('countOverlapsInBins'))
setMethod('countOverlapsInBins','GRanges',
          function(query, subject, nbins,strandspecific=FALSE) {
    if(!is(subject,'GRanges'))
        stop('subject has to be of class GRanges ..')
    if(!is.numeric(nbins))
        stop('nbins has to be of class numeric ..')
    if(!is.logical(strandspecific))
        stop('strandspecific has to be of either TRUE or FALSE ..')        
            countMat <- matrix(0, length(query), nbins)
            binsize <- floor(width(query)/nbins)
            starts <- start(query)
            if(strandspecific==TRUE){
              strands <- strand(query)
            }else{
              strands <- rep("*",length(query))
            }
            for (bin in 1:nbins) {
              startPos <- starts + binsize * (bin - 1)
              if (bin == nbins) 
                endPos <- end(query)
              else endPos <- startPos + binsize - 1
              queryBin <- GRanges(seqnames = seqnames(query), ranges = IRanges(start = startPos, 
                                                                               end = endPos),strand=strands)
              countMat[, bin] <- countOverlaps(queryBin, subject, maxgap = 0L, 
                                               type = "any")
            }
            countMat[countMat > 1] <- 1
            return(countMat)
})

