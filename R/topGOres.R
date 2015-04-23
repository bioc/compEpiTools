topGOres <- function(ids, ontology='BP', Pthr=1e-5, maxN=5000,
                     minN=5, orgdb, allEG=keys(orgdb), showTerms=NULL) {
    if(!is.character(ids))
        stop('ids has to be of class character ...')
    if(!all(ontology %in% c('BP','CC','MF')))
        stop('ontology has to be one of: CC, BP, MF ...')
    if(!is.numeric(Pthr))
        stop('Pthr has to be of class character ...')
    if(!is.numeric(maxN))
        stop('maxN has to be of class character ...')
    if(!is.numeric(minN))
        stop('minN has to be of class character ...')
    if(!is(orgdb,'OrgDb'))
        stop('orgdb has to be of class OrgDb ...')
    if(species(orgdb) != 'Homo sapiens' && species(orgdb) != 'Mus musculus' && species(orgdb) != 'Drosophila melanogaster')
        stop('only org.Hs.eg.db, org.Mm.eg.db and org.Dm.eg.db are supported ...')
    if(!is(allEG, 'character'))
        stop('allEG has to be of class character ...')
    if(!is.null(showTerms) && !is.numeric(showTerms))
      stop('showTerms has to be either NULL or of class numeric ...')
    
    EGvec <- unique(c(allEG, ids))
    EGvec <- rep(0, length(allEG))
    names(EGvec) <- allEG
    EGvec[ids] <- 1
    EGvec <- factor(EGvec)
    if(length(levels(EGvec)) == 1)
        stop('the allEG gene universe cannot be the same as the provided ids ...')
    if(species(orgdb) == 'Homo sapiens') dbname <- 'org.Hs.eg.db'
    if(species(orgdb) == 'Mus musculus') dbname <- 'org.Mm.eg.db'
    if(species(orgdb) == 'Drosophila melanogaster') dbname <- 'org.Dm.eg.db'
    
    GOdata <- new("topGOdata", ontology=ontology, allGenes=EGvec,
                  annot=annFUN.org, mapping=dbname, nodeSize=minN)
    res <- runTest(GOdata, algorithm="classic", statistic="fisher")
    resdf <- GenTable(GOdata, classic=res, topNodes=50)
    ann.genes <- genesInTerm(GOdata, whichGO=resdf$GO.ID)
    Genes <- sapply(ann.genes,function (x) {intersect(x,ids)})  
    Genes <- sapply(Genes,function(x) {c(paste(x,sep="",collapse = "/"))})
    Genes <- as.data.frame(Genes)
    resdf <- cbind(resdf,Genes)
    resdf$Term <- Term(GOTERM[resdf$GO.ID])
    pvals <- resdf$classic
    pvals <- as.numeric(sub('<', '', pvals))
    inds <- which(resdf$Annotated < maxN & pvals < Pthr)
    if(length(inds) == 0) return(NULL)
    resdf <- resdf[inds,]
    if(is.numeric(showTerms))
    {
      if(showTerms > 50)
        showTerms <- 50
      if(showTerms > nrow(resdf))
        showTerms <- nrow(resdf)
      enrich_GO <- resdf[1:showTerms,]
      enrich_GO$Term <- factor(enrich_GO$Term, levels=enrich_GO$Term)
      enrich_GO$P_val <- -log10(as.numeric(enrich_GO$classic))
      g.plot <- ggplot(enrich_GO, aes(Term, Significant, fill=P_val)) + geom_bar(stat='identity') + coord_flip()+ ylab("Number of Genes")
      return(list(data=resdf,enrichplot=g.plot))
    }
    return(resdf)
}

