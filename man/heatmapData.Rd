\name{heatmapData}
\alias{heatmapData}
\title{Based on a list of GRanges, determine various kind of counts
  before displaying a heatmap}
\description{Based on a list of GRanges, determine various kind of
  counts before displaying a heatmap}
\usage{heatmapData(grl, refgr=grl[[1]], useScore=rep(FALSE,
length(grl)), type, Nnorm=TRUE, Snorm=TRUE, txdb=NULL, nbins=5)}
\arguments{
  \item{grl}{list; a list of GRanges or paths to BAM files}
  \item{refgr}{GRanges; the reference set of genomic regions}
  \item{useScore}{logical; an optional boolean array of length equal to
    the length of grl}
  \item{type}{character; an array of length equal to the length of grl,
    with a combination of 'mcols', 'gr' or 'cov'}
  \item{Nnorm}{logical; whether to perform library size normalization,
    only applied if some of the element in type are equal to 'cov'}
  \item{Snorm}{logical; whether to perform normalization based on the
    refgr widths, only applied if some of the element in type are equal
    to 'cov'}
  \item{nbins}{numeric; the number of bins the ranges in refgr have to
    be divided into}
  \item{txdb}{an object of class TxDb}
}
\details{
  The functions is used to determine various kind of counts for each
  object in grl in each range of refgr and is typically used to prepare
  the input for the heatmapPlot method.
  
  The type of counts is determined through the corresponding type setting.
  If type is mcols, the counts are expected to be pre-calculated and
  available in the mcols of the correponding grl GRanges.
  If type is gr, the corresponding grl GRanges (gr) is considered and the
  counts are the number of occurrencies of gr for each ranges of refgr; 
  if nbins is greater than 1 and type is gr, the counts are determined
  for each bin of each range of refgr.
  A score (the lower, the more significant) can be provided in the first
  column of the mcols of gr; the minimum score over the gr ranges
  associated to every given refgr range is determined and stored in the
  corresponding column of the scoreMat output matrix.
  
  If type is cov, the corresponding grl has to be a path to a BAM file,
  and the counts are the coverage within each range of refgr; 
  if nbins is greater than 1 and type is cov, the counts are determined
  for each bin of each range of refgr.
  If Nnorm is TRUE and type is cov, the counts are divided by the
  million mapped reads in the BAM file.
  If Snorm is TRUE and type is cov, the counts are divided by the range
  width in bp.
  
  If a TxDb is provided, the presence of an intron or exon is registered
  for each range of refgr; intron is assigned 0.6, exon 0.4, and they will
  be rendered using the heatmapPlot function as red and pink, respectively.
  If nbins is greater than 1 and a TxDb is provided, the presence of an
  intron or exon is registered for each bin of each range of refgr.

  The bam files have to be associated to the corresponding index .bai
  files. Please refer to the documentation of samtools on how to create them.
}
\value{
  A list of two items, matList and scoreMat is returned.
  matList: if a TxDb is not provided, matList is a list of length equal
  to the length of grl; each item of the list is a matrix with number of
  rows equal to the number of ranges in refgr, and number of columns
  equal to nbins; if a TxDb is provided, matList is a list of length
  equal to the length of grl + 2 is returned; the two extra items
  contain the count for introns and exons.
  scoreMat: if useScore is all FALSE then scoreMat is set to NULL,
  otherwise it is a matrix whose number of rows is equal to the length
  of refgr and the number of columns is equal to the length of grl; row
  Nr and column Nc contain the minimum score of mcols(grl[[Nc]])[,1] for
  the ranges overlapping with refgr[Nr], if any (0 otherwise).
}
\references{
  http://genomics.iit.it/groups/computational-epigenomics.html
}
\examples{
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
isActiveSeq(txdb) <- c(TRUE, rep(FALSE, length(isActiveSeq(txdb)) - 1))
TSSpos <- TSS(txdb)
gr <- TSSpos[1:5]
start(gr) <- start(gr) - 1000
end(gr) <- end(gr) - 600
extgr <- GRanges(seqnames(gr), ranges=IRanges(start(gr) - 1000, end(gr) + 1000))
data <- heatmapData(grl=list(ChIPseq= gr), refgr=extgr, type='gr')
restoreSeqlevels(txdb)
}

