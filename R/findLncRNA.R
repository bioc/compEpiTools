findLncRNA <- function(k4me3gr, k4me3bam, k4me1bam, k79bam=NA,
                       k36bam=NA, RNAseqbam=NA, sizeLNC=10000,
                       extDB=NULL, txdb, org=NULL, Qthr=0.95) {
    if(!is.character(k4me3bam) || !file.exists(k4me3bam))
        stop('k4me3bam has to point to an existing BAM file ...')
    if(!is.character(k4me1bam) || !file.exists(k4me1bam))
        stop('k4me1bam has to point to an existing BAM file ...')
    if(!is.na(k79bam) && !file.exists(k79bam))
        stop('k79bam has to be either NA or a character
		pointing to an existing BAM file ...')
    if(!is.na(k36bam) && !file.exists(k36bam))
        stop('k36bam has to be either NA or a character
		pointing to an existing BAM file ...')
    if(!is.na(RNAseqbam) && !file.exists(RNAseqbam))
        stop('RNAseqbam has to be either NA or a character
		pointing to an existing BAM file ...')
    if(!is.numeric(sizeLNC))
        stop('sizeLNC has to be of class numeric ...')
    if(sizeLNC < 0)
        stop('sizeLNC has to be a positive integer ...')
    if(!is.null(extDB) && !is(extDB,'GRanges'))
	stop('extDB has to be either NULL or an object of class GRanges ...')
    if(!is(txdb,'TxDb'))
        stop('txdb has to be an object of class TxDb ...')
    if(!is.null(org) && !is(org,'BSgenome'))
        stop('org has to be either NULL or an object of class BSgenome ...')
    if((any(!is.na(c(k79bam, k36bam, RNAseqbam))) && Qthr > 0) && is.null(org))
        stop('filtering was required defining one of k79bam, k36bam,
		RNAseqbam and setting Qthr>0, but org was not set ...')
    if(!is.numeric(Qthr))
        stop('Qthr has to be an object of class numeric ...')
    if(Qthr < 0 || Qthr > 1)
        stop('Qthr has to be a numeric in [0,1] ...')

                                        # extended gb: genebodies plus 10kb on either side
    gb<- unlist(transcriptsBy(txdb, by='gene'))
    starts<- start(gb)- 1e4
    starts[starts < 1]<- 1
    start(gb)<- starts
    end(gb)<- start(gb) + 1e4 - 1

                                        # get k4me3 peaks outside extended gb
    gbinds<- which(countOverlaps(k4me3gr, gb) > 0)
    if(length(gbinds) > 0) k4me3gr<- k4me3gr[-gbinds]

                                        # get distal k4me3 peaks where k4me1 signal is lower than k4me3
    k4me1s<- GRcoverage(k4me3gr, k4me1bam, Nnorm=TRUE, Snorm=FALSE)
    k4me3s<- GRcoverage(k4me3gr, k4me3bam, Nnorm=TRUE, Snorm=FALSE)
    inds<- which(k4me3s > k4me1s)
    if(length(inds) > 0) k4me3gr<- k4me3gr[inds]
    else return(NULL)

                                        # force distal TSS-like k4me3 peaks to sizeLNC kb
                                        # downstream or upstream the peak mid point
    k4me3mp<- GRmidpoint(k4me3gr)
    lncROIleft<- k4me3mp
    start(lncROIleft)<- start(lncROIleft)- sizeLNC +1
    lncROIright<- k4me3mp
    end(lncROIright)<- end(lncROIright)+ sizeLNC -1

                                        #counting reads in ROIs
    bamfiles<- c(k4me3bam, k4me1bam, k79bam, k36bam, RNAseqbam)
    lncROIscoreLeft<- list()
    lncROIscoreRight<- lncROIscoreLeft
    for(i in 1:length(bamfiles)) {
        if(is.na(bamfiles[i])) {
            lncROIscoreLeft[[i]]<- rep(NA, length(lncROIleft))
            lncROIscoreRight[[i]]<- rep(NA, length(lncROIright))
        }
        else {
            lncROIscoreLeft[[i]]<- GRcoverage(lncROIleft, bamfiles[i],
                                              Nnorm=TRUE, Snorm= FALSE)
            lncROIscoreRight[[i]]<- GRcoverage(lncROIright, bamfiles[i],
                                               Nnorm=TRUE, Snorm= FALSE)
        }
    }
    names(lncROIscoreLeft)<- c('H3K4me3', 'H3K4me1', 'H3K79me2',
                               'H3K36me3', 'RNAseq')
    names(lncROIscoreRight)<- names(lncROIscoreLeft)

                                        # exporting
    lncRNAmat<- NULL
    for(i in 1:length(lncROIscoreLeft)) lncRNAmat<-
        cbind(lncRNAmat, signif(lncROIscoreLeft[[i]],3))
    for(i in 1:length(lncROIscoreRight)) lncRNAmat<-
        cbind(lncRNAmat, signif(lncROIscoreRight[[i]],3))
    colnames(lncRNAmat)<- c(paste(names(lncROIscoreLeft),'_up', sep=''),
                            paste0(names(lncROIscoreLeft), '_down'))
    rownames(lncRNAmat)<- GRanges2ucsc(k4me3mp)

    if(!is.null(extDB)) {
        mdbROIscore<- list()
        for(i in 1:length(bamfiles)) {
            if(is.na(bamfiles[i])) mdbROIscore[[i]]<- rep(NA, length(extDB))
            else mdbROIscore[[i]]<- GRcoverage(extDB, bamfiles[i],
                                               Nnorm=TRUE, Snorm=FALSE)
        }
        names(mdbROIscore)<- c('k4me3', 'k4me1', 'k79me2', 'k36me3', 'RNAseq')
        mdbROImat<- NULL
        for(i in 1:length(mdbROIscore)) mdbROImat<-
            cbind(mdbROImat, signif(mdbROIscore[[i]],3))
        mdbROImat<- cbind(mdbROImat, matrix(NA, nrow(mdbROImat), ncol(mdbROImat)))
        rownames(mdbROImat)<- GRanges2ucsc(extDB)
        lncRNAmat<- rbind(lncRNAmat, mdbROImat)
    }

    if(any(!is.na(c(k79bam, k36bam, RNAseqbam))) && Qthr > 0) {
                                        # profiling 100k random genomic regions each of 10kb
        chrsize<- seqlengths(seqinfo(org))
                                        # grep to exclude chrM and randoM chromosomes
        extraChr<- grep('m', names(chrsize), ignore.case=TRUE)
        if(length(extraChr) > 0) chrsize=chrsize[-extraChr]
        randomRegionsN<- round(2e5*chrsize/sum(as.numeric(chrsize)))
        Rregions<- NULL
        for(i in 1:length(randomRegionsN)) {
            Rstarts<- sample(1:chrsize[i], randomRegionsN[i])
            Rregions<- rbind(Rregions, data.frame(names(chrsize)[i],
                                                  Rstarts, Rstarts+1e4))
        }
        randomRegions100k10k<- GRanges(Rle(Rregions[,1]),
                                       ranges=IRanges(start=Rregions[,2], end=Rregions[,3]))
        refToAvoid<- c(lncROIleft, lncROIright, extDB, ignore.mcols=TRUE)
        inds<- which(countOverlaps(refToAvoid, randomRegions100k10k) > 0)
                                        # excluding random regions overlapping to previously tested regions
        if(length(inds)>0) randomRegions100k10k<- randomRegions100k10k[-inds]
                                        # keeping 1e5 random regions
        randomRegions100k10k<- randomRegions100k10k[1:1e5]
        randomRegionsList<- list()
        for(i in 1:length(bamfiles)) {
            if(is.na(bamfiles[i])) randomRegionsList[[i]]<- NA
            else randomRegionsList[[i]]<- GRcoverage(randomRegions100k10k,
                                                     bamfiles[i], Nnorm=TRUE, Snorm=FALSE)
        }
                                        # 95th thresholds
        thr95th<- sapply(randomRegionsList, quantile, Qthr, na.rm=TRUE)
        inds<- numeric()
        for(j in 3:5) inds<- c(inds, which(lncRNAmat[,j] >= thr95th[j] |
                                           lncRNAmat[,j+5] >= thr95th[j]))
        if(length(inds) > 0) lncRNAmat<- lncRNAmat[sort(unique(inds)),]
        else lncRNAmat<- NULL
    }

    return(lncRNAmat)
}
