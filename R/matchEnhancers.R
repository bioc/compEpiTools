setGeneric('matchEnhancers', function(enhGR, minD=2e4, maxD=2e5, txdb, EG2GS,
                                      upstream=2000, downstream=1000, TFGR=NULL)
           standardGeneric('matchEnhancers'))
setMethod('matchEnhancers','GRanges', function(enhGR, minD=2e4, maxD=2e5, txdb,
                                               EG2GS, upstream=2000, downstream=1000, TFGR=NULL){

    if(!is.numeric(minD))
        stop('minD has to be of class numeric ...')
    if(!is.numeric(maxD))
        stop('maxD has to be of class numeric ...')
    if(!is.numeric(upstream))
        stop('upstream has to be of class numeric ...')
    if(!is.numeric(downstream))
        stop('downstream has to be of class numeric ...')
    if(minD <= 0)
        stop('minD has to be positive ...')
    if(maxD <= 0)
        stop('minD has to be positive ...')
    if(downstream <= 0) stop('downstream has to be positive ...')
    if(maxD <= 0) stop('minD has to be positive ...')
    if(!is(txdb, "TxDb")) stop('txdb has to be of class TxDb ..')
    if(!is(EG2GS,"OrgDb")) stop('EG2GS has to be of class OrgDb ..')

    TSSgr <- TSS(txdb)
                                        # if a TF binding promoters is given, TF peaks at promoters are the reference
    if(!is.null(TFGR)) {
                                        # TF bound promoters
        mP <- GRangesInPromoters(Object=TFGR, txdb=txdb,
                                 upstream=upstream, downstream=upstream)
        mP <- distanceFromTSS(mP, txdb, EG2GS)
    }
    else mP <- TSSgr # .. otherwise, TSS are the reference
    names(mP) <- NULL

                                        # distance from closest enhancer
    mPdist2enh <- as.data.frame(distanceToNearest(x=mP, subject=enhGR))
                                        # if closest enhancer is more distant than 200Kb -> XmP
                                        #  (references not having enhancers within maxD)
    XmP.GR <- mP[which(mPdist2enh$distance > maxD)]

                                        # if closest enhancer is less than 200Kb then
                                        # check that other TSS are not in between (direct enhancer)
    closeInds <- which(mPdist2enh$distance <= maxD & mPdist2enh$distance >= minD)
    mP2closeEnh <- mPdist2enh[closeInds,]
    midpMs <- start(GRmidpoint(mP[mP2closeEnh$queryHits]))
    midpEs <- start(GRmidpoint(enhGR[mP2closeEnh$subjectHits]))
    SEdf <- cbind(midpMs, midpEs)
    inds <- which((SEdf[,2]- SEdf[,1])<0)
    SEdf[inds,] <- SEdf[inds, 2:1]
    M2Egr <- GRanges(seqnames(mP[mP2closeEnh$queryHits]),
                     ranges=IRanges(start=SEdf[,1], end=SEdf[,2]))
    mcols(M2Egr) <- mcols(mP[mP2closeEnh$queryHits])
    M2Egs.genesOv <- findOverlaps(M2Egr, TSSgr)
    mPid <- as.character(M2Egr$nearest_gene_id[queryHits(M2Egs.genesOv)])
    TSSid <- as.character(names(TSSgr)[subjectHits(M2Egs.genesOv)])
                                        # the TSS between the reference and the closer enhancer is a TSS of the same gene
    sameGene <- which(mPid == TSSid)
                                        # these won't be considered as False Positives
    M2Egs.genesOv <- M2Egs.genesOv[-sameGene]
    mP2discard <- unique(queryHits(M2Egs.genesOv))
    mP2closeEnhSel <- mP2closeEnh[-mP2discard,]

    if(is.null(TFGR)) {
        EmP.E.GR <- enhGR[mP2closeEnhSel$subjectHits]
        EmP.mP.GR <- mP[mP2closeEnhSel$queryHits]
        return(list(XP=XmP.GR, EP.E=EmP.E.GR, EP.P=EmP.mP.GR))
    }

                                        # if the closest direct enhancer is not TF bound then EmP
    EboundCount <- countOverlaps(enhGR[mP2closeEnhSel$subjectHits], TFGR)
    EmP.E.GR <- enhGR[mP2closeEnhSel$subjectHits[EboundCount == 0]]
    EmP.mP.GR <- mP[mP2closeEnhSel$queryHits[EboundCount == 0]]

                                        # if the closest direct enhancer is TF bound then mEmP
    mEmP.mE.GR <- enhGR[mP2closeEnhSel$subjectHits][EboundCount > 0]
    mEmP.mP.GR <- mP[mP2closeEnhSel$queryHits[EboundCount > 0]]

    resList <- list(XmP=XmP.GR, EmP.E=EmP.E.GR, EmP.mP=EmP.mP.GR,
                    mEmP.mE=mEmP.mE.GR, mEmP.mP=mEmP.mP.GR)
    return(resList)
})
