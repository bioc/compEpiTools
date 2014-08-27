plotStallingIndex <- function(matlist, xlimlist=NULL,
                              colors=rainbow(length(matlist))) {
    if(!is(matlist, 'list'))
        stop('matlist has to be of class list ...')
    if(!is.character(colors))
        stop('colors has to be of class characters ...')
    if(length(colors) != length(matlist))
        stop('the length of colors has to match the length of matlist ...')

    layout(matrix(1:3,1,3, byrow=FALSE))
    par('mar'=c(5,4,2,0.5))
    l <- length(matlist)
    li <- vector(length=l)
    for(i in 1:l) li[i] <- nrow(matlist[[i]])
    for(i in 1:l) matlist[[i]][matlist[[i]] == 0] <-
        min(matlist[[i]][matlist[[i]] > 0])

    if(is.null(xlimlist)){
        r <- apply(matlist[[1]],2,range)
        if(l > 1) { for(i in 2:l) r <- rbind(r,apply(matlist[[i]],2,range)) }
        xlimlist <- list(SI=c(min(r[,'SI']), max(r[,'SI'])), TSS=c(min(r[,'TSS']),
                                                                 max(r[,'TSS'])), c(min(r[,'GB']), max(r[,'GB'])))
    }
    else {
        if(length(xlimlist) != 3)
            stop('you need to supply 3 x ranges')
        if(any(sapply(xlimlist,function(x){ x[2] <= x[1] })))
            stop('the ranges supplied are invalid')
    }

    plot(sort(matlist[[1]][,'SI'], decreasing=FALSE), 100*(1:li[1])/li[1],
         main='Stalling Index', xlim=c(xlimlist$SI[1],xlimlist$SI[2]), type='l',
         lwd=2, cex.axis=1.2, xlab='', log='x', ylab='Stalling Index', col=colors[1])
    if(l > 1) {
        for(i in 2:l) { points(sort(matlist[[i]][,'SI'], decreasing=FALSE),
                               100*(1:li[i])/li[i], type='l', lwd=2, col=colors[i])
                    }
    }
    if(!is.null(names(matlist))) legend('topleft', legend=names(matlist),
                                        lwd=2, col=colors[1:l], bty='n')

    plot(sort(matlist[[1]][,'TSS'], decreasing=FALSE), 100*(1:li[1])/li[1],
         main='TSS',  xlim=c(xlimlist$TSS[1], xlimlist$TSS[2]), type='l',
         lwd=2, cex.axis=1.2, xlab='', log='x',
         ylab='Average coverage on  TSS', col=colors[1])
    if(l>1)	{
        for(i in 2:l) {
            points(sort(matlist[[i]][,'TSS'], decreasing=FALSE),
                   100*(1:li[i])/li[i], type='l', lwd=2, col=colors[i])
        }
    }

    plot(sort(matlist[[1]][,'GB'], decreasing=FALSE), 100*(1:li[1])/li[1],
         main='Gene body', xlim=c(xlimlist$GB[1],xlimlist$GB[2]), type='l',
         lwd=2, cex.axis=1.2, xlab='', log='x',
         ylab='Average coverage on  GB', col=colors[1])
    if(l>1)	{
        for(i in 2:l) {
            points(sort(matlist[[i]][,'GB'], decreasing=FALSE),
                   100*(1:li[i])/li[i], type='l', lwd=2, col=colors[i])
        }
    }
}
