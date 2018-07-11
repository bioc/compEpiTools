setGeneric('getPromoterClass', function(object, Nproc=1, org,
                                        upstream=1000, downstream=0)
  standardGeneric('getPromoterClass'))
setMethod('getPromoterClass','TxDb', function(object, Nproc=1,
                                              org, upstream=1000, downstream=0) {
  if(!is.numeric(Nproc))
    stop('Nproc has to be of class numeric ..')
  if(!is(org,"BSgenome"))
    stop('org has to be of class BSgenome.. ')
  if(!is.numeric(upstream))
    stop('upstream has to be of class numeric ..')
  if(!is.numeric(downstream))
    stop('downstream has to be of class numeric ..')
  if(sum(upstream+downstream) < 500)
    stop('sum of upstream and downstream has to be atleast 500 ..')
  
  suppressWarnings(prom <- promoters(txdb, upstream=upstream, downstream=downstream, 
                                     columns=c('tx_id','tx_name','gene_id')))
  inds <- which(start(prom)<0)
  if(length(inds)>0) start(prom)[inds]=1
  tssGR <- unique(sort(prom))
  
  # chromosome level function for parallelization
  pcChr <- function(CHR, tssGR, ORG) {
    promlen <- unique(width(tssGR))
    ucscAnnListChr <- tssGR[seqnames(tssGR) == CHR]
    CGpattern <- DNAString('CG')
    Cpattern <- DNAString('C')
    Gpattern <- DNAString('G')
    # retireving the genomic sequence
    chrseq <- ORG[[CHR]]
    promoterAvgCpGratios <- NULL
    promoterCategory <- NULL
    for(i in 1:length(ucscAnnListChr)) {
      promlen <- unique(width(tssGR))
      Pseq <- subseq(chrseq, start(ucscAnnListChr)[i], end(ucscAnnListChr)[i])
      # sliding 500bp windows each 5bp
      startPos <- seq(1, promlen - 499, 5)
      endPos <- startPos + 499
      # determining their genomic sequences
      PseqView <- Views(Pseq, startPos, endPos)
      Cs <- vcountPattern(Cpattern, PseqView)
      Gs <- vcountPattern(Gpattern, PseqView)
      PseqView_new <- DNAStringSet(PseqView)
      CpGratio <- (vcountPattern(CGpattern, PseqView_new)*500)/(Cs * Gs)
      CplusGcount <- Cs + Gs
      # the mean CpG ration is the promoter CpG ratio
      promoterAvgCpGratios[i] <- mean(CpGratio, na.rm=TRUE)
      # using the filters defined by weber et al, Nature Genet 2007
      # to define LCP (lowCG), ICP (intCG) or HCP (highCG)
      if(max(CpGratio, na.rm=TRUE) <= 0.48) {
        promoterCategory[i] <- 'lowCG'
        next
      }
      CplusGcount <- CplusGcount / 500
      HcpgInds <- which(CpGratio > 0.75)
      HcgInds <- which(CplusGcount > 0.55)
      if(length(HcpgInds) > 0 && length(intersect(HcpgInds, HcgInds)) > 0) {
        promoterCategory[i] <- 'highCG'
        next
      }
      # if not asigned to either lowCG or highCG it is assigned to intCG
      promoterCategory[i] <- 'intCG'
    }
    mcols(ucscAnnListChr)$promoterCpG <-  promoterAvgCpGratios
    mcols(ucscAnnListChr)$promoterClass <- promoterCategory
    return(ucscAnnListChr)
  }
  
  # parallelizing through the chromosomes
  chrs <- unique(as.character(seqnames(tssGR)))
  cl <- makeCluster(Nproc, type='PSOCK')
  # loading the compEpiTools and genomic sequence library on each node
  clRes <- clusterApplyLB(cl, chrs, pcChr, tssGR=tssGR, ORG=org)
  promcomb <- do.call("c",clRes)
  stopCluster(cl)
  return(promcomb)
})
