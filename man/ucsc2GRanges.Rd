\name{ucsc2GRanges}
\alias{ucsc2GRanges}
\title{Convert UCSC-formatted genomic positions into a GRanges}
\description{Convert UCSC-formatted genomic positions to a GRanges
  object}
\usage{
ucsc2GRanges(ucscPositions)
}
\arguments{
  \item{ucscPositions}{character; a vector of UCSC formatted genomic
    positions.}
}
\details{
  UCSC formatted genomic positions such as chr5:10000-35000 are
  converted into GRanges format.
}
\value{a GRanges with chr, start and end derived from the UCSC format}
\references{http://genomics.iit.it/groups/computational-epigenomics.html}
\examples{
ucsc2GRanges(c('chr1:100-500','chrX:20-1000'))
}
