\name{GRcoverage-methods}
\docType{methods}
\alias{GRcoverage}
\alias{GRcoverageInbins}
\alias{GRcoverage-methods}
\alias{GRcoverageInbins-methods}
\alias{GRcoverage,GRanges-method}
\alias{GRcoverageInbins,GRanges-method}
\title{based on a GRanges and a BAM file, returns the total coverage for
  each range, or for each bin of the range.}
\description{
  GRcoverage returns the total coverage for each range in GRanges, while
  GRcoverageInbins divides each range in equally-sized bins 
  and returns the total coverage for each bin.
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRcoverage(Object, bam, Nnorm=TRUE, Snorm=TRUE)
    
    GRcoverageInbins(Object, bam, Nnorm=TRUE, Snorm=TRUE, Nbins)
    
    where:
    
    \itemize{
      \item Object: GRanges
      \item bam: character; a path to a BAM file
      \item Nnorm: logical; whether to perform library size normalization
      \item Snorm: logical; whether to perform normalization based on
      the region width in bp
      \item Nbins: numeric; the number of equally sized bins each range
      in Object is divide into
    }
    
    For each range in the GRanges (or for each bin within the range),
    the sum of the base coverage in the range is determined.
    
    If Nnorm is TRUE, the coverage is divided by the number of aligned
    reads in the BAM file, and multiplied by 1e6.
    
    If Snorm is TRUE, the coverage is divided by the region width in bp.
    
    GRcoverage returns an array of region coverages with length equal to
    the length of the Object GRanges.
    
    GRcoverageInbins returns a matrix with nrow equal to the length of
    the Object GRanges and ncol equal to Nbins; if Nbins is equal to 1 a
    matrix of 1 column is returned.

    An increasing number of bins, does not correspond to a significant
    increase in the computation time (65sec for 43k regions and 10 bins).

    GRcoverageInbins with Nbins equal to 1 returns counts identical to
    GRcoverage but it is slightly slower (64sec vs 46sec for 43k regions).
    
    GRcoverageInbins will return NAs and a warning for the rows
    corresponding to those ranges whose width is lower than the number
    of bins.
    
    The coverage for sequences not available in the BAM file is set to 0.
    
    The BAM file has to be associated to the corresponding index .bai
    file. Please refer to the documentation of samtools on how to create
    it.
}}
\keyword{methods}
\examples{
bampath <- system.file("extdata", "ex1.bam", package="Rsamtools")
gr <- GRanges(seqnames= Rle(c('seq1','seq2')), 
	ranges=IRanges(start=c(1000, 100), end=c(2000, 1000)))
GRcoverage(Object=gr, bam=bampath, Nnorm=TRUE, Snorm=TRUE)
GRcoverageInbins(Object=gr, bam=bampath, Nnorm=TRUE, Snorm=TRUE, Nbins=4)
}
