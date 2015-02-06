setGeneric('GRcoverage', function(Object, bam, Nnorm=TRUE, Snorm=FALSE)
           standardGeneric('GRcoverage'))
setMethod('GRcoverage','GRanges', function(Object, bam,
                                           Nnorm=TRUE, Snorm=FALSE) {
    if(!is.character(bam))
        stop('bam has to be a file path of class character ..')
    if(!is.logical(Nnorm) || !is.logical(Snorm))
        stop('Nnorm and Snorm have to be of class logical ..')

                                        # checking for sequences not represented in the BAM file
    BAMseqs <- names(scanBamHeader(bam)[[1]]$targets)
    matchingSeqs <- which(as.character(seqnames(Object)) %in% BAMseqs)

    if(length(matchingSeqs) == 0) return(rep(0, length(Object)))
    param <- ApplyPileupsParam(which=Object[matchingSeqs], what='seq')

    coverage <- unlist(applyPileups(PileupFiles(bam),
                                    FUN=function(x) sum(x$seq), param=param))
                                        # normalizing per Million mapped reads
    if(Nnorm) {
        param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))
        nreads <- (countBam(bam, param= param)$records) / 1e6
        coverage <- coverage / nreads
    }
                                        # normalizing per bp
    if(Snorm) coverage <- coverage / width(Object[matchingSeqs])

    if(length(matchingSeqs) < length(Object)) {
                                        # coverage of 0 is set for ranges on sequences not in the BAM file
        coverageTot <- rep(0, length(Object))
                                        # coverage is only computed for sequences available in the BAM file
        coverageTot[matchingSeqs] <- coverage
    }
    else coverageTot <- coverage

    return(coverageTot)
})

