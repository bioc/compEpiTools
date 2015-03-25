TSS <- function(txdb) {
    if(!is(txdb,'TxDb'))
        stop('txdb has to be an object of class TxDb ...')

    TSSr <- promoters(txdb, upstream=0, downstream=0)
    ann <- select(txdb, keys=transcripts(txdb)$tx_name,
                  keytype='TXNAME', columns=c('TXNAME','GENEID'))
    ind <- which(duplicated(TSSr$tx_name)==TRUE)
    if(length(ind) >0)
    {
      TSSr <- TSSr[-c(ind)]
      ann <- ann[-c(ind),]
    }
    rownames(ann) <- ann[,1]
    names(TSSr) <- ann[TSSr$tx_name,'GENEID']
                                        # take min pos for - and max for + strand genes
    ind <- which(as.character(strand(TSSr))=='+')
    pos_start <- start(TSSr)[ind]
    end(TSSr)[ind] <- pos_start

    ind <- which(as.character(strand(TSSr))=='-')
    pos_end <- end(TSSr)[ind]
    start(TSSr)[ind] <- pos_end
    return(TSSr)
}
