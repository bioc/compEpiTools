setGeneric('GRanges2ucsc', function(Object)
           standardGeneric('GRanges2ucsc'))
setMethod('GRanges2ucsc','GRanges', function(Object) {

    ucscf <- paste0(as.character(seqnames(Object)), ':',
                    start(Object), '-', end(Object))

    return(ucscf)
})
