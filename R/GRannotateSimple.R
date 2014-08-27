setGeneric('GRannotateSimple', function(gr, txdb, upstream=2000, downstream=1000, plot=TRUE)
           standardGeneric('GRannotateSimple'))
setMethod('GRannotateSimple','GRanges', function(gr, txdb, upstream=2000,
                                                 downstream=1000, plot=TRUE) {
    if(!is(txdb,'TxDb'))
        stop('txdb has to be an object of class TxDb ...')
    if(!is.numeric(upstream))
        stop('upstream has to be an object of class numeric ...')
    if(!is.numeric(downstream))
        stop('downstream has to be an object of class numeric ...')
    if(!is.logical(plot))
      stop('plot has to be either TRUE or FALSE ...')
    

    grList <- list()

                                        # promoter
    suppressWarnings(promoterGR <- promoters(txdb, upstream=upstream,
                                             downstream=downstream))
    promoterCounts <- countOverlaps(gr, promoterGR)
    promoterInds <- which(promoterCounts > 0) # gr ranges in promoters
    if(length(promoterInds) > 0) {
        grList$promoter <- gr[promoterInds]
        gr <- gr[-promoterInds]
    }
    else grList$promoter <- NA

                                        # intragenic
    trGR <- transcripts(txdb)
    intragenicCounts <- countOverlaps(gr, trGR)
                                        # gr ranges in genebodies AND not in promoters
    intragenicInds <- which(intragenicCounts > 0)
    if(length(intragenicInds) > 0) {
        grList$intragenic <- gr[intragenicInds]
        gr <- gr[-intragenicInds]
    }
    else grList$intragenic <- NA

                                        # intergenic
    if(length(gr) > 0) {
        grList$intergenic <- gr
    }
    else grList$intergenic <- NA
    
    if(plot)
    {
      anno_type <- sapply(grList, length)
      labels <- paste0(names(anno_type),"\n", anno_type)
      pie(anno_type, main="Genomic Annotation", labels=labels)
    }

    return(grList)
})
