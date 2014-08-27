stallingIndex <- function(BAMlist,inputList=NULL, peakGRlist, peakGB=FALSE,
                          genesList, transcriptDB, countMode="transcript", upstream=300,
                          downstream=300, cutoff=600, elongationOffset=0) {

    if(!is(BAMlist, 'list'))
        stop('BAMlist has to be an object of class list ...')
    if(!all(sapply(BAMlist, file.exists)))
        stop('BAMlist has to contain path(s) to existing BAM file(s) ...')
    if(!is.null(inputList) && !is(inputList, 'list'))
        stop('inputList has to be either NULL or an object of class list ...')
    if(!is.null(inputList) && !all(sapply(inputList, file.exists)))
        stop('inputList has to be either NULL or
      contain path(s) to existing BAM file(s) ...')
    if(!is.null(inputList) && length(BAMlist) != length(inputList))
        stop('the number of inputfiles does not match the number of BAMfiles ...')
    if(!is(peakGRlist, 'list'))
        stop('peakGRlist has to be an object of class list ...')
    if(!all(sapply(peakGRlist, function(x) is(x, 'GRanges'))))
        stop('peakGRlist has to contain GRanges object(s) ...')
    if(!is.logical(peakGB))
        stop('peakGB has to be an object of class logical ...')
    if(length(BAMlist) != length(peakGRlist))
        stop('the number of items in the peakGRlist does not
    match the number of BAMfiles ...')
    if(!is(genesList, 'list')) stop('genesList has to be
    an object of class list ...')
    if(!all(sapply(genesList, function(x) is(x, 'character'))))
        stop('genesList has to contain characters ...')
    if(length(BAMlist) != length(genesList))
        stop('the number of items in the genesList does not match
      the number of BAMfiles ...')
    if(!is(transcriptDB,"TxDb")) stop('transcriptDB has
    to be of class TxDb ..')
    if(!(countMode %in% c("transcript","average","gene")))
        stop('Counting mode must be one of the following:
      gene, transcript or average ...')
    if(!is.numeric(upstream))
        stop('upstream has to be of class numeric ..')
    if(!is.numeric(downstream))
        stop('downstream has to be of class numeric ..')
    if(!is.numeric(cutoff))
        stop('cutoff has to be of class numeric ..')
    if(!is.numeric(elongationOffset))
        stop('elongationOffset has to be of class numeric ..')


    getBound <- function(peakGR, tssGR, up, down, GB=FALSE, gbGR=NULL) {
        print(strand(tssGR)[1:5])
        print(strand(tssGR)[1:5] == '+')
        iplus <- which(as.character(strand(tssGR)) == '+')
        start(tssGR[iplus]) <- start(tssGR[iplus]) - up
        end(tssGR[iplus]) <- end(tssGR[iplus] + down)

        start(tssGR[-iplus]) <- start(tssGR[-iplus]) - down
        end(tssGR[-iplus]) <- end(tssGR[-iplus] + up)

        ovTSS <- findOverlaps(peakGR, tssGR)
        indsTSS <- subjectHits(ovTSS)

        if(GB) {
            ovGB <- findOverlaps(peakGR, gbGR)
            indsGB <- subjectHits(ovGB)
            inds <- intersect(indsTSS, indsGB)
        }
        else inds <- indsTSS

        return(inds)
    }

                                        # DEFINE GRANGES for PROMOTERS AND GENE BODIES; THE TSS GR IS BETWEEN
                                        # TSS-DOWNSTREAM AND TSS+DOWNSTREAM (STRAND-SPECIFIC)
    tssGR <- promoters(transcriptDB, upstream=upstream,
                       downstream=downstream, columns=c("gene_id","tx_name","tx_id"))
    gbGR <- transcripts(transcriptDB, columns=c("gene_id", "tx_name", "tx_id"))
    if(countMode == "gene"){
       tssGR$gene_id <- as.character(tssGR$gene_id)
       gbGR$gene_id <- as.character(gbGR$gene_id)
       gbplus <- gbGR[which(as.character(strand(gbGR)) == '+')]
       gbminus <- gbGR[which(as.character(strand(gbGR)) == '-')]

                                        # FIND GENES WITH MAX END AND MIN START IN - STRAND (LONGEST ISOFORMS)
       stopminus <- tapply(end(gbminus), gbminus$gene_id, max)
       startminus <- tapply(start(gbminus), gbminus$gene_id, min)
       a <- seqnames(gbminus)
       a <- as.character(a)
       names(a) <- gbminus$gene_id
       aU <- a[unique(gbminus$gene_id)]
       seqnames <- aU[names(stopminus)]
       gene_id <- names(stopminus)
                                        # DEFINE GRANGE FOR - STRAND
       minus <- GRanges(Rle(seqnames), IRanges(startminus,stopminus),
                        gene_id=gene_id, strand='-')

                                        # FIND GENES WITH MAX END AND MIN START IN + STRAND (LONGEST ISOFORMS)
       stopplus <- tapply(end(gbplus), gbplus$gene_id, max)
       startplus <- tapply(start(gbplus), gbplus$gene_id, min)
       a <- seqnames(gbplus)
       a <- as.character(a)
       names(a) <- gbplus$gene_id
       aU <- a[unique(gbplus$gene_id)]
       seqnames <- aU[names(stopplus)]
       gene_id <- names(stopplus)
                                        # DEFINE GRANGE FOR + STRAND - GENEBODY
       plus <- GRanges(Rle(seqnames), IRanges(startplus,stopplus),
                       gene_id=gene_id, strand='+')

                                        # DEFINE GLOBAL GRANGE - GENEBODY
       suppressWarnings(gbGR<-c(plus,minus))
       gbGR <- sort(gbGR)

                                        # DO THE SAME ON TRANSCRIPT
       tssplus <- tssGR[which(as.character(strand(tssGR)) == '+')]
       tssminus <- tssGR[which(as.character(strand(tssGR)) == '-')]

                                        # FIND THE RIGHTMOST TSS ON THE - STRAND
       valuesminus <- tapply(end(tssminus), tssminus$gene_id, max)
       a <- seqnames(tssminus)
       a <- as.character(a)
                                        # END OF TSS IS UPSTREAM+DOWNSTREAM TO THE LEFT
       stop_tss <- valuesminus - (upstream + downstream + 1)
       names(a) <- tssminus$gene_id
       aU <- a[unique(tssminus$gene_id)]
       seqnames <- aU[names(valuesminus)]
       gene_id <- names(valuesminus)
                                        # DEFINE GRANGE FOR - STRAND TSS
       minus <- GRanges(Rle(seqnames), IRanges(stop_tss,valuesminus),
                        gene_id=gene_id, strand='-')

                                        # FIND THE LEFTMOST TSS ON THE + STRAND
       valuesplus <- tapply(start(tssplus), tssplus$gene_id, min)
       a <- seqnames(tssplus)
       a <- as.character(a)
                                        # END OF TSS IS UPSTREAM+DOWNSTREAM TO THE RIGHT
       stop_tss <- valuesplus + (upstream + downstream + 1)
       names(a) <- tssplus$gene_id
       aU <- a[unique(tssplus$gene_id)]
       seqnames <- aU[names(valuesplus)]
       gene_id <- names(valuesplus)
                                        # DEFINE GRANGE FOR + STRAND TSS
       plus <- GRanges(Rle(seqnames), IRanges(valuesplus, stop_tss),
                       gene_id=gene_id, strand='+')

                                        # DEFINE GLOBAL GRANGE - TSS
       suppressWarnings(tssGR<-c(plus,minus))
       tssGR <- sort(tssGR)

       names(tssGR) <- tssGR$gene_id
       tssGR <- tssGR[gbGR$gene_id]
   }

                                        # TAKE ONLY GENES WITH GB LONGER THAN CUTOFF
    selGenes <- width(gbGR) > cutoff
    tmp <- tssGR[selGenes]
    tssGR.ref <- tssGR[selGenes]
    gbGR.ref <- gbGR[selGenes]
    tssGR.ref$gene_id <- as.character(tssGR.ref$gene_id)
                                        # REMOVE NA
    tssGR.ref <- tssGR.ref[which(!is.na(tssGR.ref$gene_id))]
    gbGR.ref$gene_id <- as.character(gbGR.ref$gene_id)
    gbGR.ref <- gbGR.ref[which(!is.na(gbGR.ref$gene_id))]

                                        # DEFINE STARTS AND STOPS FOR GB IN BOTH STRANDS: ON +:
                                        # [START + DOWNSTREAM, END]; ON -: [START, END-DOWNSTREAM]
                                        # ELONGATION OFFSET IS THE QUEUE AFTER THE END OF THE GENE
    start <- ifelse(as.character(strand(gbGR.ref)) == '+',
                    start(gbGR.ref) + downstream, start(gbGR.ref) - elongationOffset)
    stop <- ifelse(as.character(strand(gbGR.ref)) == '+',
                   end(gbGR.ref) + elongationOffset, end(gbGR.ref) - downstream)
                                        # CREATE NEW GRANGE FOR GB
    refinedGB <- GRanges(seqnames(gbGR.ref), IRanges(as.vector(start),
                                                     as.vector(stop)), strand=strand(gbGR.ref),
                         elementMetadata=mcols(gbGR.ref))
                                        # REMOVE GENES WITH SEQNAMES CONTAINING "_"
    inds1 <- grep("_", as.character(seqnames(tssGR.ref)))
    inds2 <- grep("_", as.character(seqnames(refinedGB)))
    removeind <- union(inds1,inds2)
    if(length(removeind)>0)
    {
      tssGR.ref <- tssGR.ref[-c(removeind),]
      refinedGB <- refinedGB[-c(removeind),]
    }

                                        # REMOVE DUPLICATES IF YOU ARE COUNTING GENE BY GENE
    if(countMode == "gene"){
        names(tssGR.ref) <- tssGR.ref$gene_id
       	names(refinedGB) <- refinedGB$elementMetadata.gene_id
        tssGR.ref <- tssGR.ref[unique(tssGR.ref$gene_id),]
        refinedGB <- refinedGB[unique(refinedGB$elementMetadata.gene_id),]
    }

    res <- list()
    for(i in 1:length(BAMlist)){
        if(peakGB) {
            indPeak <- getBound(peakGRlist[[i]], tssGR.ref, up=2000 - upstream,
                                down=1000 - downstream, GB=TRUE, gbGR=refinedGB)
        }
        else {
            indPeak <- getBound(peakGRlist[[i]], tssGR.ref, up=2000 - upstream,
                                down=1000 - downstream, GB=FALSE, gbGR=NULL)
        }

        tssGR.tmp <- tssGR.ref[indPeak]
        GB.tmp <- refinedGB[indPeak]

        tssGR.tmp <- tssGR.tmp[tssGR.tmp$gene_id %in% genesList[[i]]]
        GB.tmp <- GB.tmp[GB.tmp$elementMetadata.gene_id %in% genesList[[i]]]

        if(length(tssGR.tmp) > 0) {

                                        # COUNT READS AND OBTAIN STALLING INDEX (TSS/GB)
            if(!is.null(inputList)){
                tssCount <- GRcoverage(tssGR.tmp, BAMlist[[i]],
                                       Nnorm=TRUE, Snorm=TRUE)
                gbCount <- GRcoverage(GB.tmp, BAMlist[[i]], Nnorm=TRUE, Snorm=TRUE)
                tssCount_input <- GRcoverage(tssGR.tmp, inputList[[i]],
                                             Nnorm=TRUE, Snorm=TRUE)
                gbCount_input <- GRcoverage(GB.tmp, inputList[[i]],
                                            Nnorm=TRUE, Snorm=TRUE)
                tssCount <- sapply(tssCount-tssCount_input, max, 0)
                gbCount <- sapply(gbCount-gbCount_input, max, 0)
                stalling_index <- tssCount/(gbCount + 0.01)
            }
            else {
                tssCount <- GRcoverage(tssGR.tmp, BAMlist[[i]],
                                       Nnorm=TRUE, Snorm=TRUE)
                gbCount <- GRcoverage(GB.tmp, BAMlist[[i]], Nnorm=TRUE, Snorm=TRUE)
                stalling_index <- tssCount/(gbCount + 0.01)
            }
            res[[i]] <- cbind(tssCount,gbCount,stalling_index)
            colnames(res[[i]]) <- c('TSS','GB','SI')

            if(countMode == "transcript"){
                rownames(res[[i]]) <- mcols(tssGR.tmp)$tx_name
            }

                                        # FOR THIS COUNT MODE, AVERAGE OVER ISOFORMS
            if(countMode == "average"){
          	rownames(res[[i]]) <- mcols(tssGR.tmp)$gene_id
                tmpTSS <- split(res[[i]][,1], rownames(res[[i]]))
                geneTSS <- lapply(tmpTSS, mean)
                tmpGB <- split(res[[i]][,2], rownames(res[[i]]))
                geneGB <- lapply(tmpGB, mean)
                tmpSI <- split(res[[i]][,3], rownames(res[[i]]))
                geneSI <- lapply(tmpSI, mean)
                res[[i]] <- cbind(unlist(geneTSS), unlist(geneGB), unlist(geneSI))
                colnames(res[[i]]) <- c('TSS','GB','SI')
                rownames(res[[i]]) <- names(geneSI)
            }

            if(countMode == "gene"){
                rownames(res[[i]]) <- mcols(tssGR.tmp)$gene_id
            }
     	}
        else{ res[[i]] <- NULL}
    }

    if(!is.null(names(BAMlist))) names(res)=names(BAMlist)
    return(res)
}

