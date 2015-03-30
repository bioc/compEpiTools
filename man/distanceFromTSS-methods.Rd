\name{distanceFromTSS-methods}
\docType{methods}
\alias{distanceFromTSS}
\alias{distanceFromTSS-methods}
\alias{distanceFromTSS,GRanges-method}
\title{Returns the GRanges annotated with info about the closer TSS}
\description{
  For each GRanges region it decorates the GRanges with extra columns containing info about the closer TSS.
}
\section{Methods}{
  \describe{
    {
      To be used in this form:
      
      distanceFromTSS(Object, txdb, EG2GS=NULL)
      
      where:
      
      \itemize{
	\item Object: GRanges
	\item txdb: TxDb
	\item EG2GS: an object of class OrgDb; like org.Mm.eg.db, org.Hs.eg.db (use the exact name of object)
      }
      
      The method returns a GRanges with additional columns.
      If EG2GS is NULL three columns are appended containing info for the gene with the closer TSS:
      \itemize{
	\item nearest_tx_name: the transcript id
	\item distance_fromTSS: the distance in bp
	\item nearest_gene_id: gene id
      }
      
      If EG2GS is provided, gene symbols are also included as additional column.
    }
}}
\keyword{methods}
\examples{
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
isActiveSeq(txdb) <- c(TRUE, rep(FALSE, length(isActiveSeq(txdb))-1))
TSSpos <- TSS(txdb)
gr <- TSSpos[1:5]
start(gr) <- start(gr)-1000
end(gr) <- end(gr)-600
mcols(gr) <- NULL
distanceFromTSS(Object=gr, txdb=txdb, EG2GS=NULL)
restoreSeqlevels(txdb)
}
