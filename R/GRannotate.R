setGeneric('GRannotate', function(Object, txdb, EG2GS, upstream=2000,
                                  downstream=1000, userAnn=NULL) standardGeneric('GRannotate'))
setMethod('GRannotate','GRanges', function(Object, txdb, EG2GS,
                                           upstream=2000, downstream=1000, userAnn=NULL) {
    if(min(width(Object)) < 1 || max(width(Object)) > 1)
        stop('GRannotate: A GRanges with width 1 is expected ..')
    if(!is(txdb, "TxDb"))
        stop('txdb has to be of class TxDb ..')
    if(!is(EG2GS,"AnnDbBimap"))
        stop('EG2GS has to be of class AnnDbBimap ..')
    if(!is.numeric(upstream))
        stop('upstream has to be of class numeric ..')
    if(!is.numeric(downstream))
        stop('downstream has to be of class numeric ..')
    if(!is.null(userAnn) && !is(userAnn,'GRangesList'))
        stop('userAnn has to ba an object of class GRangesList ...')
    if(length(unique(names(userAnn))) !=  length(userAnn))
        stop('userAnn has to have unique names ...')

    ObjectNoMcols <- Object
    mcols(ObjectNoMcols) <- NULL
                                        # distance from closer TSS and its annotation
    TSSdist <- distanceFromTSS(Object=ObjectNoMcols, txdb=txdb, EG2GS=EG2GS)

                                        # defining gene bodies and promoters
    gb <- transcriptsBy(txdb, by='gene')
    suppressWarnings(prom <- promoters(gb, upstream=upstream, downstream=downstream))
    gb <- unlist(gb)
    prom <- unlist(prom)
    ind <- which(duplicated(gb$tx_name)==TRUE)
    if(length(ind) > 0)
    {
      gb <- gb[-c(ind)]
      prom <- prom[-c(ind)]
    }

                                        # cutting promoter regions downstream the TSS from genebody regions
    shortinds <- which(width(gb) <= downstream)
    gbNotProm <- gb[-shortinds]
    plusStr <- which(as.logical(strand(gbNotProm) == '+'))
    minusStr <- which(as.logical(strand(gbNotProm) == '-'))
    start(gbNotProm[plusStr]) <- start(gbNotProm[plusStr]) + downstream
    end(gbNotProm[minusStr]) <- end(gbNotProm[minusStr]) - downstream

    loc <- rep('intergenic', length(Object))
    locEG <- rep(NA, length(Object))
    locTX <- rep(NA, length(Object))
    locGS <- rep(NA, length(Object)
                 )
                                        # looking for the overlap of Object ranges with genebodies
    gbInGR <- findOverlaps(query=gbNotProm, subject=Object,
                           type='any', select='all')
    qHits <- queryHits(gbInGR)
    sHits <- subjectHits(gbInGR)
    EG2GS <- as.list(EG2GS)
    gbHitsEG <- tapply(names(gbNotProm[qHits]), INDEX=as.factor(sHits),
                       FUN=paste, collapse=';')
    gbHitsTX <- tapply(gbNotProm[qHits]$tx_name, INDEX=as.factor(sHits),
                       FUN=paste, collapse=';')
    gbHitsGS <- tapply(EG2GS[names(gbNotProm[qHits])], INDEX=as.factor(sHits),
                       FUN=paste, collapse=';')
    gbHitsL <- tapply(rep('genebody', length(sHits)), INDEX=as.factor(sHits),
                      FUN=paste, collapse=';')
    if(length(gbHitsEG) > 0) locEG[as.numeric(names(gbHitsEG))] <- gbHitsEG
    if(length(gbHitsTX) > 0) locTX[as.numeric(names(gbHitsTX))] <- gbHitsTX
    if(length(gbHitsGS) > 0) locGS[as.numeric(names(gbHitsGS))] <- gbHitsGS
    if(length(gbHitsL) > 0) loc[as.numeric(names(gbHitsL))] <- gbHitsL

                                        # looking for the overlap of Object ranges with promoters
    promInGR <- findOverlaps(query=prom, subject=Object,
                             type='any', select='all')
    qHits <- queryHits(promInGR)
    sHits <- subjectHits(promInGR)
    promHitsEG <- tapply(names(prom[qHits]), INDEX=as.factor(sHits),
                         FUN=paste, collapse=';')
    promHitsTX <- tapply(prom[qHits]$tx_name, INDEX=as.factor(sHits),
                         FUN=paste, collapse=';')
    promHitsGS <- tapply(EG2GS[names(prom[qHits])], INDEX=as.factor(sHits),
                         FUN=paste, collapse=';')
    promHitsL <- tapply(rep('promoter', length(sHits)), INDEX=as.factor(sHits),
			FUN=paste, collapse=';')
    if(length(promHitsEG) > 0) {
        inds <- as.numeric(names(promHitsEG))
        NAinds <- which(is.na(locEG[inds]))
        nonNAinds <- which(!is.na(locEG[inds]))
        locEG[inds[NAinds]] <- promHitsEG[NAinds]
        locEG[inds[nonNAinds]] <-
            paste(locEG[inds[nonNAinds]], promHitsEG[nonNAinds], sep=';')
    }
    if(length(promHitsTX) > 0) {
        inds <- as.numeric(names(promHitsTX))
        NAinds <- which(is.na(locTX[inds]))
        nonNAinds <- which(!is.na(locTX[inds]))
        locTX[inds[NAinds]] <- promHitsTX[NAinds]
        locTX[inds[nonNAinds]] <-
            paste(locTX[inds[nonNAinds]], promHitsTX[nonNAinds], sep=';')
    }
    if(length(promHitsGS) > 0) {
        inds <- as.numeric(names(promHitsGS))
        NAinds <- which(is.na(locGS[inds]))
        nonNAinds <- which(!is.na(locGS[inds]))
        locGS[inds[NAinds]] <- promHitsGS[NAinds]
        locGS[inds[nonNAinds]] <-
            paste(locGS[inds[nonNAinds]], promHitsGS[nonNAinds], sep=';')
    }
    if(length(promHitsL) > 0) {
        inds <- as.numeric(names(promHitsL))
        INTinds <- which(loc[inds]=='intergenic')
        nonINTinds <- which(loc[inds]!='intergenic')
        loc[inds[INTinds]] <- promHitsL[INTinds]
        loc[inds[nonINTinds]] <-
            paste(loc[inds[nonINTinds]], promHitsL[nonINTinds], sep=';')
    }

                                        # assembling results and appending to Object mcols
    res <- data.frame(location=loc, location_tx_id=locTX,
                      location_gene_id=locEG, location_gene_symbol=locGS, stringsAsFactors=FALSE)
    mcols(Object) <- data.frame(mcols(Object), mcols(TSSdist), res, stringsAsFactors=FALSE)

                                        # appending user defined annotations
    if(!is.null(userAnn)) {
        for(i in 1:length(userAnn)) {
            userOv <- rep(0, length(Object))
            annOv <- findOverlaps(Object, userAnn[[i]])
            if(length(queryHits(annOv)) > 0)
                userOv[unique(queryHits(annOv))] <- 1
            mcols(Object)[[names(userAnn)[i]]] <- userOv
        }
    }

    return(Object)
})

