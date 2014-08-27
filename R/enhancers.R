setGeneric('enhancers', function(gr, txdb, upstream=2000,
                                 downstream=1000, CGIgr=NULL) standardGeneric('enhancers'))
setMethod('enhancers','GRanges', function(gr, txdb, upstream=2000,
                                          downstream=1000, CGIgr=NULL) {
    if(!is(txdb,'TxDb'))
        stop('txdb has to be an object of class TxDb ...')
    if(!is.numeric(upstream))
        stop('upstream has to be an object of class numeric ...')
    if(!is.numeric(downstream))
        stop('downstream has to be an object of class numeric ...')
    if(!is.null(CGIgr) && !is(CGIgr,'GRanges'))
        stop('CGIgr has to be either NULL or an object of class GRanges ...')
    suppressWarnings(promoterGR<- promoters(txdb, upstream=upstream, downstream=downstream))
                                        # identification of gr GRanges not overlapping with promoters
    res<- countOverlaps(gr, promoterGR)
    gr<- gr[res == 0]
                                        # optional identification of gr GRanges not overlapoing with CpG islands
                                        # might be useful to discard possible unannotated TSS
    if(!is.null(CGIgr)) {
        res<- countOverlaps(gr, CGIgr)
        gr<- gr[res == 0]
    }

    return(gr)
})
