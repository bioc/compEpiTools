setGeneric('GRannotateSimple', function(gr, txdb, upstream=2000, downstream=1000)
           standardGeneric('GRannotateSimple'))
setMethod('GRannotateSimple','GRanges', function(gr, txdb, upstream=2000,
                                                 downstream=1000) {
    if(!is(txdb,'TxDb'))
        stop('txdb has to be an object of class TxDb ...')
    if(!is.numeric(upstream))
        stop('upstream has to be an object of class numeric ...')
    if(!is.numeric(downstream))
        stop('downstream has to be an object of class numeric ...')
  
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
    
    ind <- which(is.na(grList)=="TRUE")
    if (length(ind) > 0)
      {
        grList <- grList[-c(ind)]
        return(grList)
      }
    anno_type <- sapply(grList, length)
    cols=c("brown","#ddaa00","beige")
    labels <- paste0(names(anno_type),"\n", anno_type)
    pie(anno_type, main="Genomic Annotation", labels=labels, col=cols)
    return(grList)
})
