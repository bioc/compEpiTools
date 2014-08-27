\name{unionMaxScore-methods}
\docType{methods}
\alias{unionMaxScore}
\alias{unionMaxScore-methods}
\alias{unionMaxScore,GRanges-method}
\title{GRanges method to perform union of peaks keeping the score of the
  most significant peak}
\description{
  the union between two GRanges associated to scores is performed, and
  the most significant score for the overlapping regions is returned
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    unionMaxScore(gr1, gr2, score1=mcols(gr1)[,1], score2=mcols(gr2)[,1])
    
    where:
    
    \itemize{
      \item gr1: GRanges
      \item gr2: GRanges
      \item score1: numeric; the scores (or any kind of score) associated to gr1
      \item score2: numeric; the scores (or any kind of score) associated to gr2
    }
    
    The method returns the union GRanges as in union(gr1, g2). For each
    range, if it resulted from overlapping ranges, the maximum score for
    those ranges is returned. In case pvalues (P) are adopted, they
    should be transformed in a score with something like -log10(P)
}}
\keyword{methods}
\examples{
	gr1 <- GRanges(Rle(c('chr1','chr1')), 
		ranges=IRanges(start=c(100,200), end=c(150,250)), pvalue=c(200, 100))
	gr2 <- GRanges(Rle(c('chr1','chr1')), 
		ranges=IRanges(start=c(80,400), end=c(300,500)), pvalue=c(50, 150))
	unionMaxScore(gr1, gr2)
}
