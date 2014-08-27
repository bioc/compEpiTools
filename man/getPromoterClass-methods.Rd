\name{getPromoterClass}
\docType{methods}
\alias{getPromoterClass}
\alias{getPromoterClass-methods}
\alias{getPromoterClass,TxDb-method}
\title{Determining the CpG promoter class and the average CpG content}
\description{
  Determining the CpG promoter class (Low, Intermediate or High CpG
  Content: lowCG, intCG or highCG, respectively) and the average CpG
  content for each entry of a transcriptDB.
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    getPromoterClass(txdb, Nproc=1, org, upstream=1000, downstream=0)
    
    where:
    
    \itemize{
      \item{txdb:}{An object of class TxDb}
      \item{Nproc:}{numeric; the number of processors to be used; one chr is run for each processor}
      \item{org:}{an object of class BSgenome}
      \item{upstream:}{numeric; number of bp upstream transcription start sites defining upstream limit of promoters}
      \item{downstream}:{numeric; number of bp downstream transcription start sites defining downstream limit of promoters}
    }
    According to Weber M et al, Nature Genet 2007: the CpG content is
    determined as (W*CpG)/(C*G), where W is the window size and CpG, C and
    G are the number of CpG, C and G occurrences, respectively. Here W is
    set to 500bp. The CplusG content is (C*G)/W. The promoter CpG class is
    determined sliding a 500bp window in 1Kb upstream regions, with step of 5
    bp. If the maximum CpG content for a given promoter is < 0.48, the
    promoter is assigned a lowCG. If the maximum CpG content is > 0.75 and
    CplusG>0.55, the promoter is assigned a highCG. The remaining
    promoters are assigned a intCG. The average CpG ratio is the average
    of the CpG ratios for all the windows.
    
    A \link{GRanges} object decorated with the promoterClass and promoterCpG data is returned.
}}
\keyword{methods}
\examples{
require(BSgenome.Mmusculus.UCSC.mm9)
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
isActiveSeq(txdb) <- c(rep(FALSE,20), TRUE, rep(FALSE, 14))
allpromoter <- getPromoterClass(txdb, Nproc=1, org=Mmusculus)
restoreSeqlevels(txdb)
}

