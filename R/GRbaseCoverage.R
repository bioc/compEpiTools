setGeneric('GRbaseCoverage', function(Object, bam, Nnorm=FALSE)
           standardGeneric('GRbaseCoverage'))
setMethod('GRbaseCoverage','GRanges', function(Object, bam, Nnorm=FALSE) {
    if(!is.character(bam))
        stop('bam has to be a file path of class character ..')
    if(!is.logical(Nnorm)) stop('Nnorm has to be of class logical ..')

                                        # sequences not represented in the BAM file
    BAMseqs <- names(scanBamHeader(bam)[[1]]$targets)
    matchingSeqs <- which(as.character(seqnames(Object)) %in% BAMseqs)

    if(length(matchingSeqs) == 0) return(sapply(width(Object),
                 function(x) rep(0, x))) # return a list of 0 base-coverage

    param <- ApplyPileupsParam(which=Object[matchingSeqs], what='seq')
    coverage <- applyPileups(PileupFiles(bam), FUN=function(x) x, param=param)
    widths <- width(Object[matchingSeqs])
    covList <- list()
    starts <- start(Object[matchingSeqs])
    for(i in 1:length(Object[matchingSeqs])) {
        covx <- coverage[[i]]
        cvec <- rep(0, widths[i])
                                        # positions with coverage >0, returned by applyPileups
        inds <- covx$pos-starts[i] + 1
        cvec[inds] <- colSums(covx$seq)
        covList[[i]] <- cvec
    }

    if(Nnorm) {
        param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))
        nreads <- (countBam(bam, param= param)$records) / 1e6
        covList <- lapply(covList, function(x) x / nreads)
    }

    if(length(matchingSeqs) < length(Object)) {
                                        # base-coverage of 0 is set for ranges on sequences not in the BAM file
        coverageTot <- sapply(width(Object), function(x) list(rep(0, x)))
                                        # coverage is only computed for sequences available in the BAM file
        coverageTot[matchingSeqs] <- covList
    }
    else coverageTot <- covList

    return(coverageTot)
})

