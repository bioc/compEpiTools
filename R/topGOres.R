topGOres <- function(ids, ontology='BP', Pthr=1e-5, maxN=5000,
                     minN=5, orgdb, allEG=keys(orgdb)) {
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
    if(species(orgdb) != 'Homo sapiens' && species(orgdb) != 'Mus musculus')
        stop('only org.Mm.eg.db and org.Hs.eg.db are supported ...')
    if(!is(allEG, 'character'))
        stop('allEG has to be of class character ...')

    EGvec <- unique(c(allEG, ids))
    EGvec <- rep(0, length(allEG))
    names(EGvec) <- allEG
    EGvec[ids] <- 1
    EGvec <- factor(EGvec)
    if(length(levels(EGvec)) == 1)
        stop('the allEG gene universe cannot be the same as the provided ids ...')
    if(species(orgdb) == 'Homo sapiens') dbname <- 'org.Hs.eg.db'
    if(species(orgdb) == 'Mus musculus') dbname <- 'org.Mm.eg.db'
    GOdata <- new("topGOdata", ontology=ontology, allGenes=EGvec,
                  annot=annFUN.org, mapping=dbname, nodeSize=minN)
    res <- runTest(GOdata, algorithm="classic", statistic="fisher")
    resdf <- GenTable(GOdata, classic=res, topNodes=50)
    pvals <- resdf$classic
    pvals <- as.numeric(sub('<', '', pvals))
    inds <- which(resdf$Annotated < maxN & pvals < Pthr)
    if(length(inds) == 0) return(NULL)
    resdf <- resdf[inds,]

    return(resdf)
}

