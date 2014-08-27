setGeneric('GRmidpoint', function(Object) standardGeneric('GRmidpoint'))
setMethod('GRmidpoint', 'GRanges', function(Object) {

    mp <- round(start(Object) + end(Object)) / 2
    start(Object) <- mp
    end(Object) <- mp

    return(Object)
})

