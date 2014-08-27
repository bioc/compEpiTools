\name{GRanges2ucsc-methods}
\docType{methods}
\alias{GRanges2ucsc}
\alias{GRanges2ucsc-methods}
\alias{GRanges2ucsc,GRanges-method}
\title{A GRanges method to convert ranges information into UCSC format}
\description{
  Chr assignments, start and end positions are converted into UCSC
  format, in the form chr1:100-500
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRanges2ucsc(Object)
    
    where Object is a GRanges
}}
\keyword{methods}
\examples{
gr <- GRanges(Rle(c('chr1','chr2')), 
	ranges=IRanges(start=c(1e7, 2e7), end=c(1e7+19, 2e7+19)))
GRanges2ucsc(gr)
}
