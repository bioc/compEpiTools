\name{TSS}
\alias{TSS}
\title{based on a TxDb returns a GRanges with the TSS positions for all
  transcripts}
\description{
  based on a TxDb returns a GRanges of width 1 with the TSS positions
  for all transcripts
}
\usage{TSS(txdb)}
\arguments{
  \item{txdb}{An object of class TxDb}
}
\details{
  based on a TxDb returns a GRanges of width 1 with the TSS positions
  for all transcripts
}
\value{A GRanges}
\references{http://genomics.iit.it/groups/computational-epigenomics.html}
\examples{
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
res <- TSS(txdb)
res[1:2]
}
