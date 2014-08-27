\name{countOverlapsInBins-methods}
\docType{methods}
\alias{countOverlapsInBins}
\alias{countOverlapsInBins-methods}
\alias{countOverlapsInBins,GRanges-method}
\title{given a query and a subject GRanges returns
  a matrix of counts of subject in bins of query}
\description{
  Given a query and a subject GRanges, this function returns a
  matrix with number of rows equal to the number of regions in query and
  number of columns equal to the number of bins. For each bin of each
  region, the occurrence (1) or not (0) of subject is returned.
}
\section{Methods}{
  \describe{
    {
      The countOverlaps method can be used to determine the overlap between a query and a subject GRanges.
      
      The countOverlapsInBins method add the functionality to partition query in bins.
      
      To be used in this form:
      
      countOverlapsInBins(query, subject, nbins)
      
      where:
      
      \itemize{
	\item query: GRanges
	\item subject: GRanges
	\item nbins: numeric, the number of bins
      }
      
      It returns a matrix with number of rows equal to the number of regions in query and
      number of columns equal to the number of bins.
      
      For each bin of each region 1 is assigned of subject GRanges overlap, 0 if it does not.
    }
}}
\keyword{methods}
\examples{
gr1 <- GRanges(seqnames=Rle('chr1',2), ranges=IRanges(start=c(10,100), end=c(50, 150)))
gr2 <- GRanges(seqnames=Rle('chr1',2), ranges=IRanges(start=c(2,40), end=c(15, 70)))
countOverlapsInBins(gr1, gr2, nbins=4)
}

