setGeneric('GRcoverageSummit', function(Object, bam, bamControl=NULL)
  standardGeneric('GRcoverageSummit'))
setMethod('GRcoverageSummit','GRanges', function(Object, bam, bamControl=NULL) {
  if(!file.exists(bam))
    stop('bam has to be the path to an existing file ..')
  if(!is.null(bamControl) && !file.exists(bamControl))
    stop('bamControl has to be the path to an existing file ..')
  
  if(!is.null(bamControl)) {
    blcov <- GRbaseCoverage(Object=Object, bam=bam, Nnorm=TRUE)
    blCcov <- GRbaseCoverage(Object=Object, bam=bamControl, Nnorm=TRUE)
    for(i in 1:length(blcov)) blcov[[i]]= blcov[[i]] - blCcov[[i]]
  }
  else blcov <- GRbaseCoverage(Object=Object, bam=bam, Nnorm=FALSE)
  # index of the position(s) with the max coverage
  maxPos <- lapply(blcov, function(x) which(x==max(x)))
  starts <- start(Object)
  # reconstructing the corresponding genomic positions
  for(i in 1:length(maxPos)) maxPos[[i]] <- maxPos[[i]] + starts[i] - 1
  maxPosL <- sapply(maxPos, length)
  # if multiple summits are found, take one random
  if(max(maxPosL) > 1) 
  {
    maxPos[maxPosL > 1] <-
      sapply(maxPos[maxPosL > 1], function(x) sample(x, 1))
    maxPos <- unlist(maxPos)
  }
  else maxPos <- unlist(maxPos)
  start(Object) <- maxPos
  end(Object) <- maxPos
  
  # returns the same GRanges with start and end assigned to the summit
  return(Object)
})

