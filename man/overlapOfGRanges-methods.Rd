\name{overlapOfGRanges-methods}
\docType{methods}
\alias{overlapOfGRanges}
\alias{overlapOfGRanges-methods}
\alias{overlapOfGRanges,GRanges-method}
\title{visualization of GRanges overlap}
\description{
  given a list of GRanges, all pair-wise overlap are evaluated and the
  percentage of overlapping ranges is visualized in a heatmap 
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    overlapOfGRanges(GRlist, plot=TRUE)
    
    where:
    
    \itemize{
      \item GRlist: a list of GRanges
      \item plot: logical
    }
    
    This function uses countOverlaps to count the number of shared
    ranges for all pairs of GRanges in the input list. The result is
    returned as an heatmap, white to beige to red, corresponding to
    increased overlap. The percentage overlap is added on each heatmap
    cell. For each GRanges on the column of the heatmap (GRC) it
    represents the percentage of each GRanges in the heatmap rows which
    is overlapping with GRC.
}}
\keyword{methods}
\examples{
	starts <- seq(100, 500, length.out=5)
	gr1 <- GRanges(Rle('chr1'), 
		ranges=IRanges(start=starts, end=starts+100))
	starts <- seq(300, 700, length.out=5)
	gr2 <- GRanges(Rle('chr1'), 
		ranges=IRanges(start=starts+50, end=starts+120))
	overlapOfGRanges(GRlist=list(gr1, gr2), plot=FALSE)
}
