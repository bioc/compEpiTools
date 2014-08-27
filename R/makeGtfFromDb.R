setGeneric('makeGtfFromDb', function(object, type, filename=NULL )
           standardGeneric('makeGtfFromDb'))
setMethod('makeGtfFromDb', 'TxDb', function(object, type, filename=NULL) {

    if(type != 'gene' && type != 'tx')
        stop('type has to be either "gene" or "tx" ...')
    if(!is.null(filename) && !is.character(filename))
        stop('filename has to be either NULL or of class character ...')

                                        # use method disjointExons to get a GRangeList of disjoint exons
    message( 'Retrieving data from TranscriptDb...' )

    if(type == 'gene') {
                                        # get a GRangeList of exons belonging to each gene
        exonsDB <- exonsBy(object , 'gene')
                                        # reduce each gene ( merge overlapping exons ) and transform to GRanges
        exonsDB <- unlist(reduce(exonsDB))
                                        # feature name to be passed later as annotation
        featureName <- 'gene_id "'
    }
    else if (type == 'tx') {
                                        # get a GRangeList of exons belonging to each transcript
        exonsDB <- exonsBy(object, by='tx', use.names=TRUE)
                                        # pass from GRangesList to GRanges
        exonsDB <- unlist(exonsDB)
                                        # feature name to be passed later as annotation
        featureName <- 'tx_name "'
    }

                                        # make a dataframe from the GRangeList
    exonsDF <- data.frame(
                                        # chromosome
        chromosome=seqnames( exonsDB )
                                        # source ( name of the program generating the feature )
        , source=object$packageName
                                        # feature
        , feature='exon'
                                        # start
        , start=start( exonsDB )
                                        # end
        , end=end( exonsDB )
                                        # value ( a floating point variable )
        , value='.'
                                        # strand
        , strand=strand( exonsDB )
                                        # frame ( One of '0', '1' or '2'. '0' indicates that the first base of
                                        # the feature is the first base of a codon, '1' that the second base is
                                        # the first base of a codon, and so on.. )
        , frame='.'
                                        # attribute ( A semicolon-separated list of tag-value pairs, providing
                                        # additional information about each feature. )
        , attribute=paste( featureName , names( exonsDB ) , '"' , sep='' )
	)

                                        # sorting lines according to genomic position
    message('Sorting according to genomic position...')
    idx <- order( exonsDF$chromosome , exonsDF$start )
	exonsDF <- exonsDF[idx, ]

	if(is.null(filename)) return(exonsDF)
	else {
		# write the GTF file
            write.table(
                exonsDF
                , file  = filename
                , quote = FALSE
			, sep   = '\t'
                , row.names = FALSE
                , col.names = FALSE
                )

            message('File created.')
	}
})
