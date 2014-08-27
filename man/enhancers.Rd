\name{enhancers-methods}
\docType{methods}
\alias{enhancers}
\alias{enhancers-methods}
\alias{enhancers,GRanges-method}
\title{A GRanges method to define enhancers based on H3K4me1 peaks}
\description{
  A GRanges method to define enhancers based on H3K4me1 peaks and genome annotation
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    enhancers(gr, txdb, upstream= 2000, downstream= 1000, CGIgr=NULL)
    
    where:
    
    \itemize{
      \item gr: a GRanges, typically of H3K4me1 peaks
      \item txdb: an object of class TxDB
      \item upstream: numeric; the number of bp upstream the TSS
      \item downstream: numeric; the number of bp downstream the TSS
      \item CGIgr: GRanges; optional GRanges of CpG Islands (CGI)
      
      Enhancers are defined as distal H3K4me1 peaks not overlapping with
      CGI, to avoid unannotated transcriptional units. Distal peaks are
      those peaks not overlapping with promoters. Alternative marks or
      proteins, such as H3K27ac or mediator, could be used here in place
      of H3K4me1. For example, H3K27ac would specifically allow to identify active enhancers.
    }
}}
\keyword{methods}
\examples{
	require(TxDb.Mmusculus.UCSC.mm9.knownGene)
	txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
	# loading H3K4me1 peaks as a GRanges object
  	# built based on the BED file from the GEO GSM1234488 sample
  	# limited to chr19:3200000-4000000
  	H3K4me1GR <- system.file("extdata", "H3K4me1GR.Rda", package="compEpiTools")
  	load(H3K4me1GR)
	enhancers(H3K4me1GR, txdb)
}
