\name{GRannotateSimple-methods}
\docType{methods}
\alias{GRannotateSimple}
\alias{GRannotateSimple-methods}
\alias{GRannotateSimple,GRanges-method}
\title{a GRanges method to split a GRanges in three GRanges: promoter,
  intragenic and intergenic}
\description{
  Given a GRanges and a TxDB object, the GRanges are divided according
  to their overlap to promoters (defined based on upstream and
  downstream bp from transcription start sites, TSS) and genebody
  (intragenic), while the remaining GRanges are assigned to intergenic
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRannotateSimple(gr, txdb, upstream=2000, downstream=1000)
    
    where:
    
    \itemize{
      \item gr: a GRanges
      \item txdb: an object of class TxDB
      \item upstream: numeric; the number of bp upstream the TSS
      \item downstream: numeric; the number of bp downstream the TSS
    }
}}
\keyword{methods}
\examples{
	require(TxDb.Mmusculus.UCSC.mm9.knownGene)
	txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
	gr <- GRanges(Rle(c('chr1','chr1')), 
		ranges=IRanges(start=c(100,200), end=c(150,250)))
	GRannotateSimple(gr, txdb)
}
