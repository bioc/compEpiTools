\name{GRcoverageSummit-methods}
\docType{methods}
\alias{GRcoverageSummit}
\alias{GRcoverageSummit-methods}
\alias{GRcoverageSummit,GRanges-method}
\title{Based on a GRanges and a BAM file, returns a GRanges with the
  positions of maximum coverage within each range}
\description{
Based on a GRanges and a BAM file, returns a GRangegs with the positions
of maximum coverage within each range; it can be used to identify peak
summits in ChIPseq enriched regions.
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRcoverageSummit(Object, bam)
    
    where:
    
    \itemize{
      \item Object: GRanges
      \item bam: a BAM file path
    }
    
    The method returns a GRanges with regions of width 1 pointing to the
    position of higher coverage.
    
    If multiple maxima exist in a range, one is returned at random.
    
    The bam file has to be associated to the corresponding index .bai
    file. Please refer to the documentation of samtools on how to create
    it.
}}
\keyword{methods}
\examples{
bampath <- system.file("extdata", "ex1.bam", package="Rsamtools")
gr <- GRanges(seqnames=Rle(c('seq1','seq2')), 
	ranges=IRanges(start=c(1000, 100), end=c(2000, 1000)))
GRcoverageSummit(Object=gr, bam=bampath)
}
