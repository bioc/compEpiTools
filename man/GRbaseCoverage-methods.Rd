\name{GRbaseCoverage-methods}
\docType{methods}
\alias{GRbaseCoverage}
\alias{GRbaseCoverage-methods}
\alias{GRbaseCoverage,GRanges-method}
\title{Based on a GRanges and a BAM file, returns a list of base
  coverage vectors for each range}
\description{
  Based on a GRanges and a BAM file, returns a list of base coverage
  vectors for each range
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRbaseCoverage(Object, bam, Nnorm=FALSE)
    
    where:
    
    \itemize{
      \item Object: GRanges
      \item bam: path to BAM file
      \item Nnorm: logical; whether to apply library size normalization
    }
    
    The method determines for each base in each region of the GRanges
    the number of reads in the BAM file. If Nnorm is TRUE the coverage
    is divided by million of mapped reads contained in the BAM file.
    The method returns a list of length equal to the length of the
    Object GRanges. Each list item is a vector of lenth equal to the
    length of the corresponding range width.
    The vector reports the (normalized) base coverage.
    The bam file has to be associated to the corresponding index .bai
    file. Please refer to the documentation of samtools on how to create it.
}}
\keyword{methods}
\examples{
bampath <- system.file("extdata", "ex1.bam", package="Rsamtools")
gr <- GRanges(seqnames=Rle(c('seq1','seq2')), 
	ranges=IRanges(start=c(1000, 100), end=c(2000, 1000)))
res <- GRbaseCoverage(Object=gr, bam=bampath)
str(res)
}
