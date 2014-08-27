\name{palette2d}
\alias{palette2d}
\title{build a two dimensional color palette}
\description{build a two dimensional color palette}
\usage{
palette2d(n, k, col1='white', col2='orange', col3='red') 
}
\arguments{
  \item{n}{numeric}
  \item{k}{numeric}
  \item{col1}{character; a color, the lower limit of the color palette}
  \item{col2}{character; a color, the mid point of the color palette}
  \item{col3}{character; a color, the upper limit of the color palette}
}
\details{
  A bidimensional color palette is built. First a set1 of round(n/2)
  colors ranging from col1 and col2 is defined. 
  Then a set2 of round(n/2) colors ranging from col2 and col3 is
  defined. Set1 and set2 are concatenated, as the first row of a k X n
  matrix M. Finally, for each column i of M, k colors ranging from white
  to M[1,i] are defined.

  It is used in the GRheatmapPlot function to  have a color palette
  simultaneously representing the intensity of an event (n) and its
  significance (k). See the example.
}
\value{A matrix of color codes, as described in details}
\references{http://genomics.iit.it/groups/computational-epigenomics.html}
\examples{
res <- palette2d(50, 10)
plot(rep(1:ncol(res),each=nrow(res)), rep(1:nrow(res),times=ncol(res)), col=res, 
	pch=20, xlab='intensity', ylab='significance')
}
