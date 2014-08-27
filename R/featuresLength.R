setGeneric('featuresLength', function(object, type)
           standardGeneric( 'featuresLength'))
setMethod('featuresLength', 'TxDb', function(object, type) {

                                        # controls on arguments
    if(type != 'gene' && type != 'tx')
        stop('featuresLength: type has to be either "gene" or "tx" ...')

    if(type == 'gene') {
                                        # load data from TxDb
        exonsDB<- exonsBy(object, 'gene')
        exonsDB<- unlist(reduce( exonsDB))
                                        # compute gene length
        flen<- tapply(width(exonsDB),
                      as.factor(unlist(names(exonsDB))), sum)
    }

    if(type == 'tx') {
                                        # load data from TxDb
        exonsDB<- exonsBy(object, by='tx', use.names=TRUE)
                                        # compute gene length
        flen<- sapply(width(exonsDB), sum)
    }

    return(flen)
})
