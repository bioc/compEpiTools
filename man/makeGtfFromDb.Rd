\name{makeGtfFromDb-methods}
\docType{methods}
\alias{makeGtfFromDb}
\alias{featuresLength}
\alias{makeGtfFromDb-methods}
\alias{featuresLength-methods}
\alias{makeGtfFromDb,TxDb-method}
\alias{featuresLength,TxDb-method}
\title{Utilities to transform a TxDb into a GTF file}
\description{
  Given a TxDB, the method 'makeGtfFromDb' creates and writes to file a
  GTF file of exons either at the gene level or at the transcript
  level. This methods can be used to generate a GTF file that is compliant
  with Bioconductor TxDB libraries, to be used for example with external
  tools for counting NGS reads over genes or transcripts. The method
  'featuresLength' determines the length of all the features associated to
  every gene or transctipt and returns the sum. All the exons associated
  to either a gene or a transcript are non-overlapping. This can be useful
  to determine read counts into summaries such as expression RPKMs.
}
\section{Methods}{
  \describe{
    {
      The 'makeGtfFromDb' has to be used in this form:
      
      makeGtfFromDb(object, type, filename)
      
      where:
      
      \itemize{
	\item object: TxDb
	\item type: character, either 'gene' or 'tx'
	\item filename: either NULL (a data.frame is returned) or a
	character containing the path to the file that has to be written
      }
      
      The 'featuresLength' has to be used in this form:
      
      featuresLength(object, type)
      
      where:
      
      \itemize{
	\item object: TxDb
	\item type: character, either 'gene' or 'tx'
      }
      
      It returns a vector of integers with number of items equal to the
      number genes or transcripts annotated in the TxDb given containing
      the widths of the features.
      
    }
}}
\keyword{methods}
\examples{
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
isActiveSeq(txdb) <- c(TRUE, rep(FALSE, length(isActiveSeq(txdb)) - 1))
chr1geneLengths <- featuresLength(txdb, 'gene')
res <- makeGtfFromDb(txdb, 'gene')
isActiveSeq(txdb) <- rep(TRUE, length(isActiveSeq(txdb)))
}

