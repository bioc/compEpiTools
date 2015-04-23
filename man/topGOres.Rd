\name{topGOres}
\alias{topGOres}
\title{determines GeneOntology (GO) enriched terms for a set of Entrez
  gene ids}
\description{
  determines GeneOntology (GO) enriched terms for a set of Entrez gene
  ids
}
\usage{topGOres(ids, ontology='BP', Pthr=1e-05, maxN=5000, minN=5,
orgdb, allEG=keys(orgdb),showTerms=NULL)}
\arguments{ 
  \item{ids}{character; an array of human or mouse EntrezGene ids}
  \item{ontology}{one of the three GO ontologies: BP (Biological
    Processes), CC (Cellular Components) or MF (Molecular Functions)}
  \item{Pthr}{numeric [0,1]; the p-value for an enrichment to be
    considered significant}
  \item{maxN}{numeric; only GO terms with a total up to maxN genes
    annotated on the genome are considered}
  \item{minN}{numeric; only GO terms with a minimum of minN genes
    annotated on the genome are considered}
  \item{orgdb}{An object of class OrgDb; either org.Hs.eg.db, org.Mm.eg.db
    or org.Dm.eg.db}
  \item{allEG}{character; the reference universe of EntrezGene ids, by
    default all of them}
  \item{showTerms}{numeric: the number of GO terms to plot; NULL: no plotting}  
}
\details{
  Determines GeneOntology (GO) enriched terms for a set of Entrez gene
  ids. Based on the topGO Bioconductor library. Both maxN and minN refer
  to the number of genes annotated in the reference genome for a given
  GO term (indipendently from the ids that the enrichment is being
  evaluated for). This can be used for excluding GOterm very generic or
  very specific, since these would mostly be ignored in the final
  output. This would also save time in the analysis and decrease the
  severity of the multiple testing issue.  
  
}
\value{
  A matrix containing enriched GO terms ranked by significance is
  returned, with the following columns: 
  \item{GO.ID}{GO id}
  \item{Term}{text description of the GO id}
  \item{Annotated}{total number of genes annotated with the considered
    GOterm}
  \item{Significant}{number of genes in ids for the specific GOterms}
  \item{Expected}{expected number of GOterms genes in ids in case of
    random enrichment}
  \item{Classic}{pvalue for the enrichment as reported from the topGO
    package}
  \item{Genes}{Gene ids of significantly annotated genes for each specific GOterms}
}
\references{http://genomics.iit.it/groups/computational-epigenomics.html}
\seealso{topGOdata in the topGO Bioconductor package}
\examples{
require(org.Mm.eg.db)
egs <- keys(org.Mm.egACCNUM)[1:50]
topGOres(ids=egs, Pthr=0.006, maxN=5000, minN=5, 
	orgdb=org.Mm.eg.db, allEG=keys(org.Mm.eg.db)[1:5000])
}
