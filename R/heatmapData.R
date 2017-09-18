heatmapData <- function(grl, refgr=grl[[1]], useScore=rep(FALSE, length(grl)),
                        type, Nnorm=TRUE, Snorm=TRUE, txdb=NULL, nbins=5) {
    if(!is(grl,'list') && !is(grl,'GElist'))
        stop('grl has to be of class list or GElist ...')
    if(is(grl,'list')) {
        if(any(!(sapply(grl, class) %in% c("GRanges", "character", 'GEcollection'))))
            stop('grl has to be a list of either GRanges objects,
            BAM file paths or GEcollection objects...')
    }
    if(!is(refgr,'GRanges'))
        stop('refgr has to be of class GRanges ...')
    if(length(useScore) != length(grl))
        stop('useScore has to have the same length of grl ...')
    if(!is.logical(useScore))
        stop('useScore has to be of class logical ...')
    if(!is.character(type))
        stop('type has to be of class character ...')
    if(length(type) != length(grl))
        stop('type has to have the same length of grl .. ')
    if(!is.logical(Nnorm))
        stop('Nnorm as to be of class logical ..')
    if(!is.logical(Snorm))
        stop('Sizenorm as to be of class logical ..')
    if(!is.null(txdb) && !is(txdb,'TxDb'))
        stop('txdb has to be of class TxDb ...')
    if(!is.numeric(nbins))
        stop('nbins has to be of class numeric ...')
    mcolsInds <- which(type == 'mcols')
    grInds <- which(type == 'gr')
    covInds <- which(type == 'cov')
    if(length(c(mcolsInds, grInds)) > 0 &&
       !all(sapply(grl[c(mcolsInds, grInds)], class) == 'GRanges'))
        stop('grl items with type equal to mcols or gr have
        to be of class GRanges ...')
    if(length(mcolsInds) > 0 &&
       !all(sapply(grl[mcolsInds], function(x) length(mcols(x))) == nbins))
        stop('grl items with type equal to mcols have to have
        nbins mcols columns ...')
    if(length(covInds) > 0 && !all(file.exists(unlist(grl[covInds]))))
        stop('grl items with type equal to cov have to be existing BAM files ...')

    if(is(grl,'GElist'))
        bins <- as.numeric(ncol(grl[[1]]))
    else
        {
            if(is(grl[[1]], 'GEcollection'))
                bins <- as.numeric(ncol(grl[[1]]))
            else
                bins <- nbins
        }
                                        # list to be filled with length(grl) matrices
                                        # of length(refgr) rows and nbins columns
    reslist <- list()
    types <- unique(type)

    for(i in 1:length(types)) {
        cti <- types[i]
        if(any(!cti %in% c('density', 'C','mC', 'rC','mcols','gr','cov')))
            stop('type has to be an array containing
           one of: density, C, mC, rC, mcols, gr or cov')
    }

    if(any(useScore) & is.null(txdb)) scoreMat <-
        matrix(0, length(refgr), length(grl))
    if(any(useScore) & !is.null(txdb)) scoreMat <-
        matrix(0, length(refgr), length(grl) + 2)
    else scoreMat <- NULL

    for(obji in 1:length(grl)) {
        print(names(grl)[obji])
        Gobj <- grl[[obji]]
        if(type[obji] == 'mcols') {
            data <- as.matrix(mcols(Gobj)) }
        if(type[obji] =='gr') {
            data <- countOverlapsInBins(query=refgr,
                                        subject=grl[[obji]], nbins=nbins)
                                        # a score is provided with the GRanges and is mapped to refgr
            if(useScore[obji] && ncol(mcols(grl[[obji]])) > 0) {
                sres <- findOverlaps(query=refgr, subject=grl[[obji]])
                qInds <- queryHits(sres)
                sInds <- subjectHits(sres)
                                        # the minimum score is considered
                scoreMat[unique(qInds),obji] <-
                    tapply(mcols(grl[[obji]])[sInds,1], INDEX=qInds, FUN=min)
            }
        }
        if(type[obji]== 'cov') data <- GRcoverageInbins(Object= refgr, bam= Gobj,
                   Nnorm=Nnorm, Snorm=Snorm, Nbins=bins)
        if(length(grep('C', type[obji], perl=TRUE)) > 0) {
            if(length(grep('^C', type[obji])) > 0) data <- binC(Gobj)
            if(length(grep('mC', type[obji])) > 0) data <- binmC(Gobj)
            if(length(grep('rC', type[obji])) > 0) data <- binrC(Gobj)
        }
        if(type[obji] == 'density')
            data <- binscore(Gobj)
        reslist[[obji]] <- data
    }

    names(reslist) <- names(grl)
    rownames(reslist[[1]]) <- names(refgr)
    Colsep <- cumsum(sapply(reslist, ncol))

    if(!is.null(txdb)) {
        gb <- transcripts(txdb)
        gbPlus <- gb[which(as.logical(strand(gb) == '+'))]
        gbMinus <- gb[which(as.logical(strand(gb) == '-'))]
        exons <- exons(txdb)
        exonPlus <- exons[which(as.logical(strand(exons) == '+'))]
        exonMinus <- exons[which(as.logical(strand(exons) == '-'))]
        gbPlus <- countOverlapsInBins(query=refgr, subject=gbPlus, nbins=bins)
        gbMinus <- countOverlapsInBins(query=refgr, subject=gbMinus, nbins=bins)
        exPlus <- countOverlapsInBins(query=refgr, subject=exonPlus, nbins=bins)
        exMinus <- countOverlapsInBins(query=refgr, subject=exonMinus, nbins=bins)

                                        # introns will eventually result as 0.6 and exons as 1
        trPlus <- gbPlus * 0.6 + exPlus * 0.4
        trMinus <- gbMinus * 0.6 + exMinus * 0.4

        reslist[['genes +']] <- trPlus
        reslist[['genes -']] <- trMinus
    }
    return(list(matList=reslist, scoreMat=scoreMat))
}
