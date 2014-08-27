setGeneric('unionMaxScore', function(gr1, gr2,
                                     score1=mcols(gr1)[,1], score2=mcols(gr2)[,1])
           standardGeneric('unionMaxScore'))
setMethod('unionMaxScore','GRanges', function(gr1, gr2,
                                              score1=mcols(gr1)[,1], score2=mcols(gr2)[,1]) {

    if(!is(gr2, 'GRanges'))
        stop('gr2 has to be of class GRanges ...')
    if(!is.numeric(score1))
        stop('score1 has to be of class numeric ...')
    if(!is.numeric(score2))
        stop('score2 has to be of class numeric ...')

    grU <- union(gr1, gr2)
    pvals <- matrix(NA, length(grU), 2)
                                        # determining max score over gr1
    ov1 <- findOverlaps(query= grU, subject=gr1, maxgap=0L,
                        minoverlap=1L, type='any', select='all')
    ov1 <- tapply(score1[subjectHits(ov1)],
                  INDEX=as.factor(queryHits(ov1)), FUN=max)
    pvals[as.numeric(names(ov1)), 1]= ov1
                                        # determining max score over gr2
    ov2 <- findOverlaps(query=grU, subject=gr2, maxgap=0L,
                        minoverlap=1L, type='any', select='all')
    ov2 <- tapply(score2[subjectHits(ov2)],
                  INDEX= as.factor(queryHits(ov2)), FUN=max)
    pvals[as.numeric(names(ov2)), 2]= ov2
                                        # determine max score over gr1 and gr2
    pvals <- apply(pvals, 1, function(x) max(x, na.rm=TRUE))
    mcols(grU) <- data.frame(score=pvals)

    return(grU)
})
