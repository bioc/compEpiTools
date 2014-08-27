setGeneric('GRenrichment', function(Object, bam, bamRef)
           standardGeneric('GRenrichment'))
setMethod('GRenrichment','GRanges', function(Object, bam, bamRef) {
    if(!is.character(bam))
        stop('bam has to be a file path of class character ..')
    if(!is.character(bamRef))
        stop('bamRef has to be a file path of class character ..')
    if(!file.exists(bam)) stop('bam is not an existing file ..')
    if(!file.exists(bamRef)) stop('bamRef is not an existing file ..')

                                        # checking for sequences not represented in the BAM file
    BAMseqs <- names(scanBamHeader(bam)[[1]]$targets)
    matchingSeqs <- which(as.character(seqnames(Object)) %in% BAMseqs)

    if(length(matchingSeqs) == 0) return(rep(NA, length(Object)))
    paramP <- ApplyPileupsParam(which=Object[matchingSeqs], what='seq')
    paramS <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))

                                        # bam coverage
    coverage <- unlist(applyPileups(PileupFiles(bam),
                                    FUN=function(x) sum(x$seq), param=paramP))
    nreads <- (countBam(bam, param=paramS)$records) / 1e6
    coverage <- coverage / nreads
    if(length(matchingSeqs) < length(Object)) {
                                        # coverage of NA is set for ranges on sequences not in the BAM file
        coverageTot <- rep(NA, length(Object))
                                        # coverage is only computed for sequences available in the BAM file
        coverageTot[matchingSeqs] <- coverage
    }
    else coverageTot <- coverage

                                        # bamRef coverage
    coverageRef <- unlist(applyPileups(PileupFiles(bamRef),
                                       FUN=function(x) sum(x$seq), param=paramP))
    nreadsRef <- (countBam(bamRef, param=paramS)$records) / 1e6
    coverageRef <- coverageRef / nreadsRef
    if(length(matchingSeqs) < length(Object)) {
                                        # coverage of NA is set for ranges on sequences not in the BAM file
        coverageTotRef <- rep(NA, length(Object))
                                        # coverage is only computed for sequences available in the BAM file
        coverageTotRef[matchingSeqs] <- coverageRef
    }
    else coverageTotRef <- coverageRef

                                        # computing enrichment
    enrichment <- log2(coverageTot - coverageTotRef)
                                        # if either the coverage of bam or bamref is NA -> the enrichment is NA
                                        # if the coverage of bam or bamref is identical -> the enrichment is -Inf

    return(enrichment)
})
