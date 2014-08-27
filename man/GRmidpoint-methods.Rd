\name{GRmidpoint-methods}
\docType{methods}
\alias{GRmidpoint}
\alias{GRmidpoint-methods}
\alias{GRmidpoint,GRanges-method}
\title{Returns a GRanges containing the mid point of a GRanges}
\description{
  Returns a GRanges containing the mid point of a GRanges
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRmidpoint(Object)
    
    where Object is a GRanges.
    
    A GRanges with width 1 containing the mid points of each range in
    Object is returned.
}}
\keyword{methods}
\examples{
gr <- GRanges(seqnames=Rle('chr1',2), 
	ranges=IRanges(start=c(10,100), end=c(50, 150)))
GRmidpoint(gr)
}
