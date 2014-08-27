setGeneric('distanceFromTSS', function(Object, txdb, EG2GS=NULL)
           standardGeneric('distanceFromTSS'))
setMethod('distanceFromTSS','GRanges', function(Object, txdb,
                                                EG2GS=NULL) {
    if(!is(txdb,"TxDb")) stop('txdb has to be of class TxDb ..')
    if(!is.null(EG2GS) && !is(EG2GS,"AnnDbBimap"))
        stop('EG2GS has to be either NULL or an object of class AnnDbBimap ..')
    EG2GS<- as.list(EG2GS)
    TSSpos<- TSS(txdb)
    nearestInd<- nearest(Object, TSSpos)
    nonNAinds<- which(!is.na(nearestInd))
    nearestId<- rep(NA, length(Object))
    nearestId[nonNAinds]<- TSSpos[nearestInd[nonNAinds]]$tx_name
    nearestEG<- rep(NA, length(Object))
    nearestEG[nonNAinds]<- names(TSSpos[nearestInd[nonNAinds]])
    if(!is.null(EG2GS)) nearestGS<- as.character(EG2GS[nearestEG])
    nearestDist<- rep(NA, length(Object))
    suppressWarnings( nearestDist[nonNAinds]<-
                     distance(Object[nonNAinds], TSSpos[nearestInd[nonNAinds]]) )
    if(is.null(EG2GS)) mcols(Object)<- data.frame(
        mcols(Object),
        nearest_tx_name=nearestId,
        distance_fromTSS=nearestDist,
        nearest_gene_id=nearestEG)
    else mcols(Object)<- data.frame(
        mcols(Object),
        nearest_tx_name=nearestId,
        distance_fromTSS=nearestDist,
        nearest_gene_id=nearestEG,
        nearest_gene_symbol=nearestGS, stringsAsFactors=FALSE)

    return(Object)
})

