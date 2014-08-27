\name{GR2fasta-methods}
\docType{methods}
\alias{GR2fasta}
\alias{GR2fasta-methods}
\alias{GR2fasta,GRanges-method}
\title{A GRanges method to extract and write to the disk a fasta file
  containing genomic sequences for the GRanges regions in a genome}
\description{
  Given a GRanges and a reference genome returns the sequences for all the
  ranges as in fasta format.
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GR2fasta(GR, org, fastaFile=NULL)
    
    where:
    
    \itemize{
      \item GR: GRanges
      \item org: an object of class BSgenome
      \item fastaFile: character or NULL; an optional file name for the
      file to be written on disk
    }
    For each range in the gr GRanges , the unmasked reference sequenced
    is retrieved. All the sequences are returned in fasta format,
    named by the genomic ranges in UCSC format, and optionally written on disk.
}}
\keyword{methods}
\examples{
	require(BSgenome.Mmusculus.UCSC.mm9)
	gr <- GRanges(Rle(c('chr1','chr2')), 
		ranges=IRanges(start=c(1e7, 2e7), end=c(1e7+19, 2e7+19)))
	show(GR2fasta(GR=gr, org=Mmusculus, fastaFile=NULL))
}
