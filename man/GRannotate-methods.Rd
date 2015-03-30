\name{GRannotate-methods}
\docType{methods}
\alias{GRannotate}
\alias{GRannotate-methods}
\alias{GRannotate,GRanges-method}
\title{Based on a GRanges and a TxDb, returns the GRanges with a series
  of annotations}
\description{
  Based on a GRanges and a TxDb, returns the GRanges with a series of annotations
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRannotate(Object, txdb, EG2GS, upstream=2000, downstream=1000, userAnn=NULL)
    
    where:
    
    \itemize{
      \item Object: a GRanges object with ranges of width 1
      \item txdb: TxDb
      \item EG2GS: an object of class OrgDb; like org.Mm.eg.db, org.Hs.eg.db (use the exact name of object)
      \item upstream: numeric, bp upstream the TSS to define promoters
      \item downstream: numeric, bp downstream the TSS to define promoters
      \item userAnn: either NULL or a named GRangesList containing user
      defined annotations as GRanges objects
    }
    
    The method returns a GRanges with extra columns containing the following annotations:
    \itemize{
      \item nearest_tx_name: transcript id for the gene with the closer TSS
      \item distance_fromTSS: distance in bp from the closer TSS
      \item nearest_gene_id: gene id for the gene with the closer TSS
      \item nearest_gene_symbol: gene symbol for the gene with the closer TSS
      \item location: 'intergenic' or a combination of 'promoter' and 'genebody'
      \item location_tx_id: transcript id(s) corresponding to the
      matched location, location is not 'intergenic'
      \item location_gene_id: gene id(s) corresponding to the matched
      location, location is not 'intergenic'
      \item location_gene_symbol: gene symbol(s) corresponding to the matched
      location, location is not 'intergenic'
    }
    
    This method only works with GRanges containing ranges of width 1. In
    case of annotation of GRanges of different width, such as ChIP-seq
    peaks, they should be summarized to GRanges of width 1. This could
    be easily done using the GRmidpoint method, pointing to the peaks
    mid points, or using GRcoverageSummit, pointing to the positions of
    maximum coverage provided that the BAM file is available. The
    obtained annotation only refers to these GRanges of width 1, and
    should not be referred as the annotation of a wider region. This is
    necessary to decrease the complexity of the output. Indeed, even
    given a simple GRanges of width 1, one could have multiple
    associated genomic features associated to each range. For example, a
    range could be associated at the same time to a promoter of a
    transcript and to the intron of another isoform for the same gene.
    
    Additional columns will be reported in the output GRanges if a userAnn
    GRangesList is provided (if userAnn is not NULL). For each of these
    extra-columns, 1 (or 0) are used to indicate overlap (or lack of) for
    each given Object range with at least one range of the correspondent
    userAnn GRanges. To this purpose, annotation tracks (tables) that are
    available in UCSC can be used as userAnn GRangesList. Please refer to
    the documentation of ucscTableQuery in the rtracklayer package to
    download these tables in R and convert them into a GRangesList
    object. See the example for a case in which mm9 CpG Islands are
    retrieved and provided to GRannotate to match each genomic region of
    interest (the mid points of gr) to this annotation source (CGIgr).
    Promoter regions are defined based on upstream and downstream bp from TSS.
}}
\keyword{methods}
\examples{
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
isActiveSeq(txdb) <- c(TRUE, rep(FALSE, length(isActiveSeq(txdb))-1))
require(org.Mm.eg.db)
require(rtracklayer)
TSSpos <- TSS(txdb)
gr <- TSSpos[1:5]
start(gr) <- start(gr) - 1000
end(gr) <- end(gr) - 600
mcols(gr) <- NULL
res <- GRannotate(Object=GRmidpoint(gr), txdb=txdb, EG2GS=org.Mm.eg.db, 
	upstream=2000, downstream=1000)
isActiveSeq(txdb) <- rep(TRUE, length(isActiveSeq(txdb)))
## alternatively, CGI can be incoirporated as follow:
## retrieving CGI mm9 islands from UCSC annotation tables
# session <- browserSession()
# genome(session) <- 'mm9'
# query <- ucscTableQuery(session, 'cpgIslandExt')
# CGIgr <- as(track(query), 'GRanges')
# res <- GRannotate(Object=GRmidpoint(gr), txdb=txdb, EG2GS=org.Mm.eg.db, 
#	upstream=2000, downstream=1000, userAnn=GRangesList(CGI=CGIgr))
}
