setGeneric('distanceFromTSS', function(Object, txdb, EG2GS=NULL)
           standardGeneric('distanceFromTSS'))
setMethod('distanceFromTSS','GRanges', function(Object, txdb,
                                                EG2GS=NULL) {
    if(!is(txdb,"TxDb")) stop('txdb has to be of class TxDb ..')
    if(!is.null(EG2GS) && !is(EG2GS,"OrgDb"))
        stop('EG2GS has to be either NULL or an object of class OrgDb ..')
    TSSpos<- TSS(txdb)
    nearestInd<- nearest(Object, TSSpos)
    nonNAinds<- which(!is.na(nearestInd))
    nearestId<- rep(NA, length(Object))
    nearestId[nonNAinds]<- TSSpos[nearestInd[nonNAinds]]$tx_name
    nearestEG<- rep(NA, length(Object))
    nearestEG[nonNAinds]<- names(TSSpos[nearestInd[nonNAinds]])
    if(!is.null(EG2GS))
    {
      anno_obj <- deparse(substitute(EG2GS))
      type <- metadata(txdb)[8,2]
      if (type == "Entrez Gene ID") {
        kt <- "ENTREZID"
      } else if (type =="Ensembl gene ID" || type == "Ensembl Gene ID") {
        kt <- "ENSEMBL"
      } else {
        warnings("geneID type in TranscriptDb is not supported...\t Gene information cannot be mapped...\n")
      }
      if(kt == "ENTREZID")
      {
        anno_obj_name <- sub("eg.db","egSYMBOL",anno_obj)
        annoDb <- eval(parse(text=anno_obj_name))
        EG2GS <- as.list(annoDb)
      }
      else if (kt == "ENSEMBL")
      {
        anno_obj_name <- sub("eg.db","egSYMBOL",anno_obj)
        annoDb <- eval(parse(text=anno_obj_name))
        anno_ens_name <- sub("eg.db","egENSEMBL",anno_obj)
        annoDb_ens <- eval(parse(text=anno_ens_name))
        EG2GS <- as.list(annoDb)
        EG2GS_ens <- as.list(annoDb_ens)
        names(EG2GS) <- EG2GS_ens
      }
      nearestGS<- as.character(EG2GS[nearestEG])
    }
      
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

