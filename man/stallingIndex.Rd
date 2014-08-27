\name{stallingIndex}
\alias{stallingIndex}
\title{returns a list with average read count on TSS, gene body, and
  stalling index for a number of samples}
\description{
  based on a TxDb, and a list of Pol2 ChIP-seq samples with a GRange
  corresponding to the peak call, returns a list whose elements contain
  the stalling index and average read count on TSS and gene body for a
  user-defined list of genes
}
\usage{stallingIndex(BAMlist,inputList=NULL, peakGRlist, peakGB=FALSE, 
  genesList, transcriptDB, countMode="transcript", upstream=300, 
  downstream=300, cutoff=600, elongationOffset=0)}
\arguments{
  \item{BAMlist}{List of characters; paths to BAM files containing ChIP-seq aligned reads}
  \item{inputList}{List of characters (optional); paths to BAM files
    containing ChIP-seq aligned reads of inputs relative to the ChIP-seq
    samples in BAMlist}
  \item{peakGRlist}{List of GRanges; Pol2 peaks relative to the ChIP-seq
    samples in BAMlist. Only genes with a peak on the TSS will be used
    for counting}
  \item{peakGB}{logical; whether or not to consider only genes which
    have also a Pol2 peak on the gene body. Deafult=FALSE}
  \item{genesList}{List of characters; Entrez gene IDs defining the
    genes where the stalling index will be computed}
  \item{transcriptDB}{An object of class transcriptDb}
  \item{countMode}{either 'transcript' or 'average' or 'gene'. Genomic
    units used to count reads: 'transcript' considers all transcripts
    separately, 'average' averages over different isoforms of a gene,
    'gene' considers the longest isoform only. Default='transcript'}
  \item{upstream}{numeric; upstream end of the TSS interval where reads
    will be counted. Default=300}
  \item{downstream}{numeric; downstream end of the TSS interval where
    reads will be counted. Default=300}
  \item{cutoff}{numeric; minimum length required for a gene to be
    considered for the counts. Default=600}
  \item{elongationOffset}{numeric; downstream end of the gene body
    interval where reads will be counted. Default=0}
}
\details{
  Given a set of Pol2 ChIP-seq samples, computes the average base coverage
  on TSS and gene body and their ratio (Pol2 stalling index) computed over
  specific sets of genes for a list of samples.  
  
  For each sample, reads will be counted only for genes contained in the
  corresponding element of the genesList, and only for those genes having
  a Pol2 peak (given in peakGRlist as list of GRanges) on the TSS region
  (defined through upstream and downstream). If peakGB is set to TRUE, a
  Pol2 peak on the body of a gene will be required. The lists of genes
  peakGRlist must be provided in terms of entrezIDs.
  
  Reads can be counted on individual gene isoforms
  (countMode='transcript'), averaging over isoforms
  (countMode='average') or taking the longest isoform for each gene
  (countMode='gene').
  
  The read counts will be performed on the intervals [TSS-upstream,
  TSS+downstream] (TSS) and [TSS+downstream, TES+elongationOffset]
  (genebody), where upstream, downstream, elongationOffset are
  user-defined parameters. Moreover, only genes having a length bigger
  than the user-defined parameter cutoff will be considered. 
  
  BAM files have to be associated to the corresponding index .bai
  files. Please refer to the documentation of samtools on how to create
  them.
}
\value{A list of matrices, each containing TSS counts, gene body counts and stalling index for each sample}
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
  restoreSeqlevels(txdb)
}
