setGeneric('GRsetwidth', function(gr, newWidth) standardGeneric('GRsetwidth'))
setMethod('GRsetwidth', 'GRanges', function(gr, newWidth) {
    if(!is.numeric(newWidth))
        stop('newWidth has to be of class numeric ...')

    mp <- GRmidpoint(gr)
    newS <- start(mp) - round(newWidth) / 2
                                        # width is not guarantee to be exactly newWidth
    newS[newS<0] <- 1
                                        # this might go beyond the chromosome length!
    newE <- start(mp)+ round(newWidth) / 2 - 1
    start(gr) <- newS
    end(gr) <- newE
    return(gr)
})
