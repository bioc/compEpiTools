\name{GRangesInPromoters-methods}
\docType{methods}
\alias{GRangesInPromoters}
\alias{GRangesInPromoters-methods}
\alias{GRangesInPromoters,GRanges-method}
\title{Based on a GRanges and a TxDb, subsets the GRanges to those
  regions overalpping with promoters}
\description{
Subset a GRanges returning only those ranges which are overlapping (at
least 1bp) with promoters
}
\section{Methods}{
\describe{
  This method subsets a GRanges returning only those ranges which are
  overlapping, or not, with promoters defined based on a TxDb, upstream
  and downstream distance from TSS.
  
  To be used in this form:
  
  GRangesInPromoters(Object, txdb, upstream=2000, downstream=1000, invert=FALSE)
  
  where:

  \itemize{
    \item Object: GRanges
    \item txdb: TxDb
    \item upstream: numeric
    \item downstream: numeric
    \item invert: logical; see below
  }
  
  Promoters are defined based on upstream and downstream distances from
  Transcription Start Sites (TSS). If invert if FALSE, the subset of
  Object overlapping with promoter regions is returned, if any,
  otherwise NULL is returned. If invert if TRUE, the subset of Object
  which is not overlapping with promoter regions is returned, if any,
  otherwise NULL is returned.
}}
\keyword{methods}
\examples{
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb<- TxDb.Mmusculus.UCSC.mm9.knownGene
TSSpos <- TSS(txdb)
gr <- TSSpos[1:5]
start(gr) <- start(gr) - 1000
end(gr) <- end(gr) - 600
GRangesInPromoters(Object=gr, txdb=txdb, upstream=2000, downstream=1000)
}
