\name{compEpiTools-package}
\alias{compEpiTools-package}
\alias{compEpiTools}
\docType{package}
\title{Tools for computational epigenomics}
\description{Tools for computational epigenomics}
\details{
  \tabular{ll}{
    Package: \tab compEpiTools\cr
    Type: \tab Package\cr
    Version: \tab 1.0.2 \cr
    Date: \tab 2015-02-24 \cr
    License: \tab GPL\cr
    Depends: \tab methods\cr
  }
  
  The package offers the following functionalities, divided by topic
  
  Counting reads in GRanges:
  \itemize{
    \item GRbaseCoverage: based on a GRanges and a BAM file, returns a list of base coverage vectors for each range
    \item GRcoverage: based on a GRanges and a BAM file, returns the total coverage for each range
    \item GRcoverageInbins: same as GRcoverage but dividing each range in equally-sized bins
    \item GRcoverageSummit: based on a GRanges and a BAM file, returns a GRanges with the positions of maximum coverage within each range
    \item GRenrichment: determines the enrichment over a set of genomic regions given two BAM files
    \item countOverlapsInBins: given a query and a subject GRanges returns a matrix of counts of subject in bins of query
    \item stallingIndex: computes the PolII stalling index based on number of ChIP-seq reads in promoter and genebody
  }
  
  Annotation of genomic regions:
  \itemize{
    \item TSS: based on a tXDb returns a GRanges with the TSS positions for all transcripts
    \item distanceFromTSS: based on a GRanges, returnes the GRanges annotated with info about the closer TSS
    \item GRangesInPromoters: based on a GRanges and a TxDb, subsets the GRanges to those regions overlapping with promoters
    \item GRmidpoint: returns a GRanges containing the mid point of a GRanges
    \item GRannotate: based on a GRanges and a TxDb, returns the GRanges with a series of annotations
    \item GRannotateSimple: a GRanges method to split a GRanges in three GRanges: promoter, intragenic and intergenic
    \item makeGtfFromDb: utilities to transform a TranscriptDb into a GTF file
  }
  
  Functional annotation:
  \itemize{
    \item enhancers: a GRanges method to define enhancers based on H3K4me1 peaks
    \item matchEnhancers: a GRanges method to match enhancers with putative targets sites (either TSS or TF-bound TSS)
    \item topGOres: determines GeneOntology enriched terms for a set of gene ids
    \item simplifyGOterms: simplify a list of GeneOntology terms based on the list of genes assigned to each GO term
    \item findLncRNA: identify putative long non coding RNAs (lncRNA) based on ChIP-seq chromatin features and RNAseq data
    \item getPromoterClass: determining the CpG promoter class and the average CpG content
  }
  
  Visualization:
  \itemize{
    \item heatmapData: Based on a list of GRanges, determine various kind of counts before displaying a heatmap
    \item palette2d: build a two dimensional color palette
    \item heatmapPlot: displays the heatmap based on the data from GRheatmapData
    \item plotStallingIndex: plot the PolII stalling index
  }
  
  Other:
  \itemize{
    \item GR2fasta: A GRanges method to extract and write to the disk a fasta file containing genomic sequences for the GRanges regions in a genome    
    \item overlapOfGRanges: given a list of GRanges, all pair-wise overlap are evaluated and the percentage of overlapping ranges is visualized in a heatmap
    \item GRsetwidth: set the width of a GRanges based on the mid point of each region
    \item unionMaxScore: GRanges method to perform union of peaks keeping the pvalue of the most significant peak
    \item GRanges2ucsc: a GRanges method to convert ranges information into UCSC format
    \item ucsc2GRanges: convert UCSC-formatted genomic positions into a GRanges
  }
}
\author{
  Computational Epigenomics Unit at the Center for Genomic Sciences of IIT@SEMM, Milan, Italy
  \url{http://genomics.iit.it/groups/computational-epigenomics.html}
  \email{mattia.pelizzola@gmail.com}
}

\keyword{ package }

