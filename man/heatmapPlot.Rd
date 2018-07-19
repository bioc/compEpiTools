\name{heatmapPlot}
\alias{heatmapPlot}
\title{displays the heatmap based on the data from heatmapData}
\description{
  displays heatmap of counts for a list of GRanges, typically computed
  based on the heatmapData function
}
\usage{heatmapPlot(matList, sigMat=NULL, qnorm=NULL, tnorm=NULL,
rowLab=FALSE, colLab=TRUE, margins=NULL, colors=NULL,
clusterInds=1:length(matList), dendrogram=TRUE)}
\arguments{
  \item{matList}{list; a list of matrices, all with the same number of
    rows and columns, typically the output of heatmapData}
  \item{sigMat}{matrix; a matrix of p-values with nrow equal to the nrow
    of the matList matrices and ncol equal to the length of matList}
  \item{qnorm}{array; an array with length equal to the length of
    matList containing either NULL or thresholds for quantile
    normalization, see details}
  \item{tnorm}{array; an aray with length equal to the length of matList
    containing either NULL or threshold for data normalization, see
    details}
  \item{rowLab}{logical; whether to add row labels to the heatmap, taken
    from the rownames of matList[[1]]}
  \item{colLab}{logical; whether to add column names to the heatmap,
    taken from the names of matList}
  \item{colors}{either NULL or character with an array of valid colors;
    it only works if sigMat is NULL}
  \item{margins}{either NULL or a numeric array of length 2} 
  \item{clusterInds}{either NULL or numerics, defining with matList
    items have to be used to drive the clustering of heatmap rows}
  \item{dendrogram}{logical; whether to display the dendrogram folling
    the clustering of rows}
}
\details{
  Each matrix in matList is either ranging from 0 to 1 or will be forced
  to by dividing to its maxumum.
  Alternatively, matrix normalization can be obtained using qnorm or tnorm.
  Setting qrnorm to X [0,1] for a given matList matrix will force the
  maximum of the matrix to be quantile(matrix, X).
  Setting trnorm to X for a given matList matrix will force the maximum of
  the matrix to be X.
  Using either qnorm or tnorm matrix will be finally normalized to 1 by
  dividing it by its maximum.
  
  If sigMat is not NULL it is expected to contain a list of pvalues or
  scores [0,1] for each range for each matList dataset.
  The colorscale of the heatmap will be adjusted to display the
  significance, by lightening the observation colors as a function of 
  its significance, see the example. A minimum pvalue of 1e-10 is
  forced. 50 levels of intensity (white to orange to red palette, as
  displayed in the colorscale) and 10 levels of significance (white for
  the less significance, full color according to the intensity palette, as
  displayed in the colorsclae) are considered. Plese ignore the values
  reported in the colorscale.
  If sigMat is NULL, the normalized intensity in each matList item is
  reported as it is on a white to beige to red default palette, or based
  on the colors defined in the colors argument.
  
  If margins is set to NULL the row and columns margins used to display
  labels are computed automatically, otherwise a numeric array of length
  two can be set to define them (in lines).
  
  clusterInds can be used to define which matList items drive the
  clustering of heatmap rows. If clusterInds is NULL no clustering is
  performed and no dendrogram is displayed. If clusterInds is an array
  of index in 1:length(matList), only those matList items will be used
  to determine the clustering and the dendrogram, while all matList data
  will be displayed.

  If a TxDB was provided to heatmapData beofre calling heatmapPlot, the
  last two tracks are about the overlap with exons and introns in the
  forward and reverse strand, repectively. If the default white to red
  color palette is used, and sigMat is NULL exons will be plotted in red
  and introns in pink. Rather, if sigMat is defined, introns will be in
  orange.
}
\value{
  A list with two items
  \item{data}{the normalized matrix used for the final heatmap visualization}
  \item{heatRes}{the list invisibly returned by the heatmap.2 function,
    see its documentation}
  
  Be careful that the rowInds contained in heatRes poiting to the new
  order of clustered data rows, is intended to list the reordered rows
  starting from the bottom of the heatmap.
}
\references{
  http://genomics.iit.it/groups/computational-epigenomics.html
}
\seealso{
  See Also as \code{\link{heatmap.2}}
}
\examples{
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
isActiveSeq(txdb) <- c(TRUE, rep(FALSE, length(isActiveSeq(txdb)) - 1))
TSSpos <- TSS(txdb)
names(TSSpos) <- TSSpos$tx_name
gr <- TSSpos[1:50]
start(gr) <- start(gr) - 1000
end(gr) <- end(gr) - 600
pvalues <- c(runif(20,1e-20,1e-8), runif(15,1e-4,1e-2), runif(15,0.5,1))
mcols(gr) <- pvalues
extgr <- GRanges(seqnames(gr), ranges= IRanges(start(gr) - 1000, end(gr) + 1000))
data <- heatmapData(grl=list(ChIPseq=gr), refgr=extgr, type='gr', useScore=TRUE,
  Nnorm=TRUE, Snorm=TRUE, nbins=6, txdb=txdb)
rownames(data[[1]][[1]]) <- paste(1:50, signif(pvalues,1), sep=' # ')
heatmapPlot(matList=data$matList, qnorm=NULL, tnorm=NULL, 
  rowLab=TRUE, colLab=TRUE, clusterInds=1:3)
dev.new()
heatmapPlot(matList=data$matList, sigMat=data$scoreMat, qnorm=NULL, tnorm=NULL, 
  rowLab=TRUE, colLab=TRUE, clusterInds=1:3)
restoreSeqlevels(txdb)
}
