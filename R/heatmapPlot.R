function (matList, sigMat = NULL, qnorm = NULL, tnorm = NULL, 
    rowLab = FALSE, colLab = TRUE, margins = NULL, colors = NULL, 
    clusterInds = 1:length(matList), dendrogram = TRUE) 
{
    if (!is(matList, "list")) 
        stop("matList has to be of class list ...")
    if (!all(sapply(matList, class) == "matrix")) 
        stop("matList has to be a list of matrices ...")
    if (!all(sapply(matList, ncol) == ncol(matList[[1]]))) 
        stop("matList has to be a list of matrices with\n        the same number of columns ...")
    if (any(sapply(matList, min, na.rm = TRUE) < 0)) 
        stop("matList cannot contain negative counts ...")
    if (!is.null(sigMat) && !is(sigMat, "matrix")) 
        stop("sigMat has to be of class matrix ...")
    if (!is.null(sigMat) && (nrow(sigMat) != nrow(matList[[1]]) || 
        ncol(sigMat) != length(matList))) 
        stop("sigMat does not have correct dimensions ...")
    if (!is.null(sigMat) && (min(sigMat) < 0 || max(sigMat) > 
        1)) 
        stop("sigMat has to range in [0,1] ...")
    if (!is.null(qnorm) && (!is.numeric(qnorm) || length(qnorm) != 
        length(matList))) 
        stop("qnorm has to be either NULL or an array of numeric having\n        the same length as matList ...")
    if (!is.null(tnorm) && (!is.numeric(tnorm) || length(tnorm) != 
        length(matList))) 
        stop("tnorm has to be either NULL or of an array of numeric having\n            the same length as matList ...")
    if (length(qnorm) > 1 && length(tnorm) > 1 && length(intersect(which(is.null(qnorm)), 
        which(is.null(tnorm)))) > 0) 
        warning("some matList items will be normalized according to\n        both qnorm and tnorm ...")
    li <- length(matList)
    if (names(matList)[li] == "genes -" && (!is.null(qnorm[[li]]) || 
        !is.null(tnorm[[li]]))) 
        warning("'genes +' and 'genes -' might be normalized ...")
    if (!is.logical(rowLab) || !is.logical(colLab)) 
        stop("colLab and rowLab have to be of class logical ...")
    if (!is.null(margins) && !is.numeric(margins)) 
        stop("margins has to be either NULL or a numeric array of length 2...")
    if (!is.null(colors) && !is.character(colors)) 
        stop("colors has to be either NULL or of class character...")
    if (is.numeric(margins) && length(margins) != 2) 
        stop("margins has to be either NULL or a numeric array of length 2 ...")
    if (!is.null(clusterInds) && !is.numeric(clusterInds)) 
        stop("clusterInds has to be either NULL or of class numeric ...")
    if (!is.null(clusterInds) && !all(clusterInds %in% 1:length(matList))) 
        stop("clusterInds has to be either NULL or an array containing\n        indexes to matList items ...")
    if (!is.logical(dendrogram)) 
        stop("dendrogram has to be of class logical ...")
    if (names(matList)[li] == "genes -" && !is.null(sigMat)) 
        sigMat[, (ncol(sigMat) - 1):ncol(sigMat)] <- 0
    if (!is.null(qnorm)) {
        for (i in 1:length(qnorm)) {
            if (is.null(qnorm[i])) 
                next
            maxQ <- quantile(matList[[i]], qnorm[[i]], na.rm = TRUE)
            matList[[i]][matList[[i]] > maxQ] <- maxQ
            matList[[i]] <- matList[[i]]/maxQ
        }
    }
    if (!is.null(tnorm)) {
        for (i in 1:length(tnorm)) {
            if (is.null(tnorm[i])) 
                next
            matList[[i]][matList[[i]] > tnorm[[i]]] <- tnorm[[i]]
            matList[[i]] <- matList[[i]]/tnorm[[i]]
        }
    }
    notNorm <- which(is.null(qnorm) & is.null(tnorm))
    if (length(notNorm) > 0) {
        for (i in 1:length(matList)) {
            if (max(matList[[i]], na.rm = TRUE) > 1) {
                matList[[i]] <- matList[[i]]/max(matList[[i]], 
                  na.rm = TRUE)
                matList[[i]][matList[[i]] > 1] <- 1
            }
        }
    }
    mat <- NULL
    for (i in 1:length(matList)) mat <- cbind(mat, matList[[i]])



    if (!is.null(clusterInds)) {
        clmat <- NULL
        for (i in clusterInds) clmat <- cbind(clmat, matList[[i]])
        NAcounts <- apply(mat, 1, function(x) length(which(is.na(x))))
        NAinds <- which(NAcounts > ncol(mat) * 0.3)
        if (length(NAinds) == nrow(mat)) 
            return(NULL)
        if (length(NAinds) > 0) {
            mat <- mat[-NAinds, ]
            clmat <- clmat[-NAinds, ]
        }
		clRes <- hclust(dist(clmat), method = "average")
        rowv <- as.dendrogram(clRes)
        if (dendrogram) 
            dendr <- "row"
        else dendr <- "none"
    }
    else {
        rowv <- FALSE
        dendr <- "none"
    }
    if (rowLab) {
        rowlab <- rownames(mat)
        rmar <- 8
    }
    else {
        rowlab <- rowLab
        rmar <- 4
    }
    if (colLab) {
        midCol <- round(ncol(matList[[1]])/2)
        labCol <- cumsum(sapply(matList, ncol)) - midCol
        collab <- rep("", ncol(mat))
        collab[labCol] <- names(matList)
        cmar <- 8
    }
    else collab <- colLab
    cmar <- 4
    if (is.null(margins)) 
        margins <- c(cmar, rmar)
    if (!is.null(sigMat)) {
        colmat <- palette2d(n = 30, k = 10)
        cols <- as.vector(t(colmat))
        colbreaks <- 0:300
        sigMat[sigMat < 1e-10] <- 1e-10
        scaledSig <- round(-log10(sigMat))
        scaledSig[scaledSig > 0] <- scaledSig[scaledSig > 0] - 
            1
        scaledSigMat <- NULL
        nbins <- ncol(matList[[1]])
        for (i in 1:ncol(scaledSig)) {
            for (j in 1:nbins) scaledSigMat <- cbind(scaledSigMat, 
                scaledSig[, i])
        }
        mat <- round(1 + mat * 29)
        mat <- mat + 30 * scaledSigMat
    }
    else {
        drange <- range(mat, na.rm = TRUE)
        colbreaks <- seq(drange[1], drange[2], length.out = 129)
        if (is.null(colors)) {
            col12 <- colorRampPalette(c("white", "beige"))
            col23 <- colorRampPalette(c("beige", "red"))
            cols <- c(col12(64), col23(64))
        }
        else {
            col12 <- colorRampPalette(colors)
            cols <- c(col12(128))
        }
    }
    colSep <- c(0, cumsum(sapply(matList, ncol)))
    res <- heatmap.2(x = mat, breaks = colbreaks, Colv = FALSE, 
        Rowv = rowv, dendrogram = dendr, scale = "none", col = cols, 
        density.info = "none", symkey = FALSE, trace = "none", 
        margins = margins, labRow = rowlab, labCol = collab, 
        sepcolor = "black", colsep = colSep, sepwidth = 0.05, 
        na.color = "grey")
    invisible(list(data = mat, heatRes = res))
}