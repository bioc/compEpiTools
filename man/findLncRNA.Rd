\name{findLncRNA}
\alias{findLncRNA}
\title{Identify putative long non coding RNAs (lncRNA)}
\description{Identify putative long non coding RNAs (lncRNA)
  based on ChIP-seq chromatin features and RNAseq data}
\usage{
findLncRNA(k4me3gr, k4me3bam, k4me1bam, k79bam, k36bam, RNAseqbam,
sizeLNC=10000, extDB= NULL, txdb, org=NULL, Qthr=0.95)
}
\arguments{
  \item{k4me3gr}{GRanges; the set of H3K4me3 peaks from a ChIP-seq
    experiment}
  \item{k4me3bam}{character; a path to BAM file containing H3K4me3
    ChIP-seq aligned reads}
  \item{k4me1bam}{character; a path to BAM file containing H3K4me1
    ChIP-seq aligned reads}
  \item{k79bam}{either NA or a path to BAM file containing H3K79me2
    ChIP-seq aligned reads}
  \item{k36bam}{either NA or a path to BAM file containing H3K36me3
    ChIP-seq aligned reads}
  \item{RNAseqbam}{either NA or a path to BAM file containing RNA-seq
    aligned reads}
  \item{sizeLNC}{numeric; the size of the putative lncRNA}
  \item{extDB}{GRanges; a set lncRNAs to be included in the analysis}
  \item{txdb}{an object of class TxDb}
  \item{org}{either NULL or an object of class BSgenome}
  \item{Qthr}{numeric in [0,1]; the percentile of the signal in random
    genomic regions to be considered as minumum cutoff}
}
\details{
  Putative long non coding RNAs (lncRNAs) are identified based on the
  associated chromatin features and, possibly, RNAseq signal.
  Only putative lncRNAs distal from genebodies are identified.
  Briefly, H3K4me3 peaks are used as the main mark indicating
  transcriptional activity. Only H3K4me3 peaks outside genebodies +/-
  10Kb are considered. Only peaks where the signal of H3K4me1 is lower
  than H3K4me3 are kept, to discard possible enhancer sites. Regions of
  interest of length sizeLNC (default 10Kb) are considered from the mid
  point of remaining H3K4me3 peaks, either in the forward or reverse
  direction (ROIs). The rational is that distal H3K4me3 marks could
  indicate the TSS of distal lncRNAs.

  Optionally, to increase the likelihood of having identified a
  bona-fide transcriptional unit, downstream regions (ROIs on either the
  forward or reverse strand) are evaluated for the existence of
  significant transcriptional signal. This is achieved based on the
  optional data (optional tracks) to be provided as BAM files: ChIP-seq
  for H3K79me2 or H3K36me3, or RNA-seq. Density of H3K79me2, H3K36me3
  and RNAseq reads in the remaining H3K4me3 regions and in 100K random
  regions of 10Kb each is determined (non overlapping with ROIs),
  normalized by the respective library size. ROIs where the signal of
  any of these is higher than the Qthr percentile of the random regions
  (profiled for the same mark) are considered as putative lncRNAs.

  An optional GRanges containing regions to be considered in any case
  (for example based on lists of a priori known lncRNAs) can be provided
  as extDB. These regions will be evaluated as they are, and subject to
  the same filtering procedure based on H3K79me2, H3K36me3 and RNAseq
  data, if provided.
  
  Passing at least one of H3K79me2, H3K36me3 and RNAseq while setting
  Qthr to 0 correspond to profile the ROIs (and the extDB regions if
  provided) for those optional tracks and avoid filtering based on the
  signal of random regions. All bam files have to be associated to the
  corresponding index .bai files. Please refer to the documentation of
  samtools on how to create them.
}
\value{Either NULL of a data.frame where putative lncRNAs (UCSC-format
  coordinates are reported as row names) are reported on the rows and
  the columns indicate the library-size normalized reads density of the
  following marks: H3K4me3, H3K4me1, H3K79me2, H3K36me3 and RNAseq
  reads. Reads density is reported for both up- and down-stream regions
  of width sizeLNC, see details, while for extDB regions the reads
  density is reported only for the regions as they are defined.}
\references{http://genomics.iit.it/groups/computational-epigenomics.html}
\examples{
  require(TxDb.Mmusculus.UCSC.mm9.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  # loading H3K4me3 peaks as a GRanges object
  # built based on the BED file from the GEO GSM1234483 sample
  # limited to chr19:3200000-4000000
  H3K4me3GR <- system.file("extdata", "H3K4me3GR.Rda", package="compEpiTools")
  load(H3K4me3GR)
  # pointing to Pol2 BAM file (it could be used as a replacement of the K79bam or K36bam ..)
  # BAM file from the GEO GSM1234478 sample, limited to chr19:3200000-4000000
  Pol2bam <- system.file("extdata", "Pol2.bam", package="compEpiTools")
  # pointing to H3K4me3 BAM file
  # BAM file from the GEO GSM1234483 sample, limited to chr19:3200000-4000000
  H3K4me3bam <- system.file("extdata", "H3K4me3.bam", package="compEpiTools")
  # pointing to H3K4me1 BAM file
  # BAM file from the GEO GSM1234488 sample, limited to chr19:3200000-4000000
  H3K4me1bam <- system.file("extdata", "H3K4me1.bam", package="compEpiTools")
  res <- findLncRNA(k4me3gr=H3K4me3GR, k4me3bam=H3K4me3bam, k4me1bam=H3K4me1bam, 
    k79bam=Pol2bam, k36bam=NA, RNAseqbam=NA, 
    sizeLNC=10000, txdb=txdb, org=NULL, Qthr=0)
}
