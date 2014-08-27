\name{GRenrichment-methods}
\docType{methods}
\alias{GRenrichment}
\alias{GRenrichment-methods}
\alias{GRenrichment,GRanges-method}
\title{Determines the enrichment over a set of genomic regions given two
  BAM files}
\description{
  Given a GRanges, and two BAM files, this method detetermines the
  coverage of the 1st BAM file (bam) and of the 2nd BAM file (bamRef),
  both normalized by millions of reads in the library. Subsequently, the
  enrichment is computed as log2(bam - bamRef).
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    GRenrichment(Object, bam, bamRef)
    
    where:
    
    \itemize{
      \item Object: a GRanges
      \item bam: path to a BAM file
      \item bamRef: path to a BAM file, to be used as reference
      
      This is typically useful when applied to BAM files derived from ChIP-seq
      experiments, where bam and bamRef are the result of the alignment of
      ChIP and input reads, respectively.
    }
    
    An array of length equal to the length of Object is returned, containing
    the enrichment values. The enrichment of a given genomic region in
    Object is set to NA if in either bam or bamRef there are no reads, while
    it is -Inf of the same normalized coverage is found for both bam and
    bamRef (as in the example).
    
    The bam file has to be associated to the corresponding index .bai
    file. Please refer to the documentation of samtools on how to create it.
}}
\keyword{methods}
\examples{
bampath <- system.file("extdata", "ex1.bam", package="Rsamtools")
gr <- GRanges(seqnames=Rle(c('seq1','seq2')), 
	ranges=IRanges(start=c(1000, 100), end=c(2000, 1000)))
 # bam and bamRef should not be pointing to the same file in real life ..
GRenrichment(Object=gr, bam=bampath, bamRef=bampath)
}
