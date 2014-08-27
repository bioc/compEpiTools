\name{plotStallingIndex}
\alias{plotStallingIndex}
\title{Stalling Index plots}
\description{
  generates the plot from the output of the stallingIndex function. The
  plot has 3 panels: Stalling Index, TSS and gene body} 
\usage{plotStallingIndex(matlist, xlimlist=NULL,
colors=rainbow(length(matlist)))} 
\arguments{
  \item{matlist}{List of matrices; each matrix must have a TSS, GB and
    SI column, as in the output of the function stallingIndex}
  \item{xlimlist}{List of numeric vectors (optional); ranges for the x
    axis of the three plots. The list must have 3 elements named 'SI',
    'TSS' and 'GB'. Default:NULL}
  \item{colors}{array; names of the colors used for the lines in the
    plot. Default:rainbow palette}
}
\details{
  Generates a 3-panel plot for the Stalling Index data (Stalling Index,
  TSS, and gene body).
}
\value{A plot}
\references{http://genomics.iit.it/groups/computational-epigenomics.html}
\examples{
  require(TxDb.Mmusculus.UCSC.mm9.knownGene)
  require(org.Mm.eg.db)
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  isActiveSeq(txdb) <- c(rep(FALSE,18), TRUE, rep(FALSE, length(isActiveSeq(txdb))-19))
  # pointing to Pol2 BAM file
  # BAM file from the GEO GSM1234478 sample, limited to chr19:3200000-4000000
  Pol2bam <- system.file("extdata", "Pol2.bam", package="compEpiTools")
  # loading Pol2 peaks as a GRanges object
  # built based on the BED file from the GEO GSM1234478 sample
  # limited to chr19:3200000-4000000
  Pol2GR <- system.file("extdata", "Pol2GR.Rda", package="compEpiTools")
  load(Pol2GR)
  egs <- distanceFromTSS(Pol2GR, txdb=txdb)
  egs <- unique(egs$nearest_gene_id)
  SI_matrix <- stallingIndex(BAMlist=list(Pol2bam), peakGRlist=list(Pol2GR), 
    genesList=list(egs), transcriptDB=txdb, countMode='gene')
  plotStallingIndex(SI_matrix)
  restoreSeqlevels(txdb)
}
