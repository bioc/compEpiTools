palette2d <- function(n, k, col1='white', col2='orange', col3='red') {
    if(!is.numeric(n)) stop('n has to be numeric ...')
    if(!is.numeric(k)) stop('k has to be numeric ...')
    if(!is.character(col1)) stop('col1 has to be character ...')
    if(!is.character(col2)) stop('col2 has to be character ...')
    if(!is.character(col3)) stop('col3 has to be character ...')

    col12 <- colorRampPalette(c(col1,col2))
    col23 <- colorRampPalette(c(col2,col3))
    colmat <- matrix(NA, k, n)
    for(i in 1:n) {
        cols <- c(col12(n / 2), col23(n / 2))
        wtocol <- colorRampPalette(c('white',cols[i]))
			colmat[,i] <- c(wtocol(k))
    }

    return(colmat)
}

