\name{GRsetwidth-methods}
\docType{methods}
\alias{GRsetwidth}
\alias{GRsetwidth-methods}
\alias{GRsetwidth,GRanges-method}
\title{Set the width of a GRanges based on the mid point of each region}
\description{
  Given a GRanges, this method sets the width of a GRanges based on the
  mid point of each region
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRsetwidth(gr, newWidth)
    
    where:
    
    \itemize{
      \item gr: a GRanges
      \item newWidth: a positive numeric
    }
}}
\keyword{methods}
\examples{
	gr <- GRanges(Rle(c('chr1','chr1')), 
		ranges=IRanges(start=c(100,200), end=c(150,250)))
	GRsetwidth(gr, 1000)
}
