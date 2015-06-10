setGeneric('GR2fasta', function(GR, org, fastaFile=NULL)
  standardGeneric('GR2fasta'))
setMethod('GR2fasta','GRanges', function(GR, org, fastaFile=NULL){
  if(!is(org,'BSgenome'))
    stop('org has to be of class BSgenome ...')
  if(!is.null(fastaFile) && !is.character(fastaFile))
    stop('fastaFile has to be either NULL or of class character ...')
  if(!is.null(fastaFile) && file.exists(fastaFile))
    stop('fastaFile already exists! ...')
  
  chrs <- as.character(seqnames(GR))
  returnseqs <- DNAStringSet()
  seqnames <- paste(as.character(seqnames(GR)),':', start(GR),
                    '-', end(GR), sep='')
  strandtype <- unique(as.character(strand(GR))) %in% "*"
  if(strandtype)
    strand(GR) <- "+"
  Chrs <- unique(chrs)
  chrs_org <- seqlengths(org)
  ind <- match(Chrs,names(chrs_org))
  suppressWarnings(seqlengths(GR) <- chrs_org[ind])
  GR_trim <- trim(GR)
  seqs <- getSeq(org, GR_trim)
  if(!is.null(fastaFile)){
    names(seqs) <- seqnames
    writeXStringSet(seqs, filepath=fastaFile, append=TRUE,
                    format="fasta")
  }else{
    returnseqs <- c(returnseqs,seqs)
    return(returnseqs)
  }
})
