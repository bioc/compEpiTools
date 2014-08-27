\name{matchEnhancers-methods}
\docType{methods}
\alias{matchEnhancers}
\alias{matchEnhancers-methods}
\alias{matchEnhancers,GRanges-method}
\title{GRanges method to match enhancers with putative targets sites}
\description{
  GRanges method to match enhancers with putative targets sites (either
  TSS or TF-bound TSS)
}
\section{Methods}{
  \describe{
    To be used in this form:
    
    matchEnhancers(enhGR, minD=2e4, maxD=2e5, txdb, EG2GS, TFGR=NULL)
    
    where:
    
    \itemize{
      \item enhGR: a GRanges of enhancer sites, such as those provided
      by the enhancers method
      \item minD: a positive numeric; the minimum distance between an
      enhancer and a target genomic region
      \item maxD: a positive numeric; the maximum distance between an enhancer
      and a target genomic region
      \item txdb: an object of class TxDb
      \item EG2GS: mapping between EntreGene and gene symbols; like
      org.Mm.egSYMBOL
      \item TFGR: an optional GRanges collecting genomic regions bound
      from a transcription factor which also binds promoters
    }
    
    This methods relies on a previous identification of enhancer sites,
    as the one performed by the enhancers compEpiTools method based on
    H3K4me1 peaks (or alternatively H3K27ac, mediator, ..).
    Genomic regions which are evaluated as putative target sites of
    these enhancers are either transcription start sites (TSS) or TSS
    which are bound at least by a transcription factor (TF) whose
    binding sites are provided wuth TFGR.
    
    Putative target sites of the provided enhancers are here defined
    based on the maximum and minimum distance from an enhancer. In
    addition, no additional TSS (belonging to differe genes, isoforms of
    the same gene do not count) have to be in between the enhancer and
    the reference target region.

    If a set of TF bound regions is also provided, which is supposed to
    contain binding sites at the level of promoters, this method returns
    a list with 5 items:
    \itemize{
      \item XmP: a GRanges with location of TSS with no enhancers, based
      on the distance contraints
      \item EmP.E: a GRanges of TF-unbound enhancers putatively
      associated to the TF-bound promoters EmP.mP
      \item EmP.mP: a GRanges of TF-bound promoters putatively
      associated to the TF-unbound enhancer sites EmP.E
      \item mEmP.mE: a GRanges of TF-bound enhancers putatively
      associated to the TF-bound promoter EmP.mP
      \item mEmP.mP: a GRanges of TF-bound promoters putatively
      associated to the TF-bound enhancer sites EmP.E
    }
    
    otherwise this method returns a list with 3 items 
    \itemize{
      \item XP: a GRanges with location of TSS with no enhancers, based
      on the distance contraints
      \item EP.E: a GRanges of enhancers putatively associated to the
      promoters EP.P
      \item EP.P: a GRanges of promoters putatively associated to the
      enhancer sites EP.E
    }
    
}}
\keyword{methods}
\examples{
	require(org.Mm.eg.db)
	require(TxDb.Mmusculus.UCSC.mm9.knownGene)
	txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
	# loading H3K4me1 peaks as a GRanges object
  	# built based on the BED file from the GEO GSM1234488 sample
  	# limited to chr19:3200000-4000000
  	H3K4me1GR <- system.file("extdata", "H3K4me1GR.Rda", package="compEpiTools")
  	load(H3K4me1GR)
	enh <- enhancers(H3K4me1GR, txdb)
	m.enh <- matchEnhancers(enhGR=enh, minD=2e4, maxD=2e5, 
		txdb=txdb, EG2GS=org.Mm.egSYMBOL)
}
