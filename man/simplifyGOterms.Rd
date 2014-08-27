\name{simplifyGOterms}
\alias{simplifyGOterms}
\title{simplify a list of GO terms}
\description{simplify a list of GeneOntology terms based on the list of
  genes assigned to each GO term}
\usage{
simplifyGOterms(goterms, maxOverlap= 0.8, ontology, go2allEGs)
}
\arguments{
  \item{goterms}{character; a vector of GO ids}
  \item{maxOverlap}{numeric in (0,1) see details}
  \item{ontology}{character; one of BP, MF or CC}
  \item{go2allEGs}{the species specific assignment of each GO term to
    EntrezGene ids}
}
\details{
  Given a pair of parent and child GO terms, and the Entrez Gene ids
  (genes sets) which can be assigned to them, the parent GO term is
  defined as redundant if the two gene sets overlap more than
  100*maxOverlap percent of the parent genes. Given a set goterms, this
  rule is applied to discard redundant parent terms while keeping the
  corresponding children terms.   
}
\value{A subset of the (possibly unaltered) input goterms}
\references{http://genomics.iit.it/groups/computational-epigenomics.html}
\examples{
require(org.Mm.eg.db)
simplifyGOterms(goterms=c('GO:0002320','GO:0002244'), maxOverlap= 0.1, ontology='BP', go2allEGs= org.Mm.egGO2ALLEGS)
}
