setGeneric('GRangesInPromoters', function(Object, txdb, upstream=2000,
                                          downstream=1000, invert=FALSE)
           standardGeneric('GRangesInPromoters'))
setMethod('GRangesInPromoters','GRanges', function(Object, txdb,
                                                   upstream=2000, downstream=1000, invert=FALSE) {

    if(!is(txdb, "TxDb"))
        stop('txdb has to be of class TxDb ..')
    if(!is.numeric(upstream))
        stop('upstream has to be of class numeric ..')
    if(!is.numeric(downstream))
        stop('downstream has to be of class numeric ..')
    if(!is.logical(invert))
        stop('invert has to be of class logical ...')

    txbygene <- transcriptsBy(txdb, by='gene')
    suppressWarnings(pregions <- unlist(promoters(txbygene,
                                                  upstream=upstream, downstream=downstream)))
    inds <- overlapsAny(query= Object,
			subject= pregions, minoverlap=1L, type='any')
    if(invert) {
        inds_f <- which(inds==FALSE)
        if(length(inds_f) > 0)
            res <- Object[inds_f]
        else
            res <- NULL
    }
    else {
        inds_t <- which(inds==TRUE)
        if(length(inds_t) > 0) res <- Object[inds_t]
        else res <- NULL
    }
    return(res)
})

