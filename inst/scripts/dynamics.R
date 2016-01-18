
remap <- function(object, ref = 1) {
    stopifnot(inherits(object, "MSnSetList"))
    x <- msnsets(object)
    refset <- x[[ref]]
    pca1 <- prcomp(exprs(refset), scale = TRUE, center = TRUE)
    preds <- lapply(x, function(xx) {
        ans <- predict(pca1, newdata = exprs(xx))
        ans
    })
    for (i in seq_along(x)) {
        xx <- object@x[[i]]
        exprs(xx) <- preds[[i]]
        sampleNames(xx) <-
            paste0("PC", 1:ncol(xx))
        xx <- MSnbase:::logging(xx, "remapped (PCA)")
        object@x[[i]] <- xx
    }
    exprs(object@x[[ref]]) <- pca1$x
    object@log[["remap"]] <- list(type = "PCA", prcomp = pca1, ref = ref)
    if (validObject(object))
        return(object)
}


## pdist <- function(object, pcol, distfun = dist) {
##     if (is.character(pcol)) {
##         stopifnot(pcol %in% names(pData(object)))
##         pcol <- pData(object)[, pcol]
##     }
##     pcol <- factor(pcol)
##     stopifnot(length(levels(pcol)) == 2)
##     stopifnot(length(pcol) == ncol(object))
##     lv <- levels(pcol)
##     e1 <- exprs(object)[, pcol == lv[1]]
##     e2 <- exprs(object)[, pcol == lv[2]]
##     ans <- sapply(seq_len(nrow(object)),
##                   function(i) distfun(rbind(e1[i, ], e2[i, ])))
##     names(ans) <- featureNames(object)
##     ans
## }

.pdist <- function(x, y, distfun) {
    e1 <- exprs(x)
    e2 <- exprs(y)
    sapply(seq_len(nrow(x)),
           function(i) distfun(rbind(e1[i, ], e2[i, ])))
}

pdist <- function(object, distfun = dist, simplify = TRUE) {
    stopifnot(inherits(object, "MSnSetList"))
    ## keeping common feature names
    suppressMessages(object <- commonFeatureNames(object))
    if (!.sameNbCol(object))
        stop("All MSnSets must have the same number of columns.")
    if (is.null(names(object)))
        names(object) <- seq_len(length(object))
    pw <- combn(length(object), 2)
    ans <- matrix(NA, ncol = ncol(pw),
                  nrow = nrow(object[[1]]))
    colnames(ans) <- apply(pw, 2,
                           function(ii) paste(names(object)[ii],
                                              collapse = "_"))
    rownames(ans) <- featureNames(object[[1]])
    for (k in 1:ncol(pw)) {
        i <- pw[1, k]
        j <- pw[2, k]
        ans[, k] <- .pdist(object[[i]], object[[j]], distfun)
    }
    if (ncol(ans) == 1 & simplify)
        ans <- ans[, 1]
    ans
}

## msnl <- commonFeatureNames(list(n2014 = nightingale2014, n2015 = nightingale2015))
## d <- pdist(msnl)

## d <- pdist(yeastgg2, "condition")
## n <- 5
## o <- order(d, decreasing = TRUE)
## d[o[1:n]]

## yggfoi <- FeaturesOfInterest(description = "top dists", fnames = names(d[o])[1:n])

## p <- plot2D(yeastgg2[, yeastgg2$condition == "Glu"])
## highlightOnPlot(yeastgg2[, yeastgg2$condition == "Glu"],
##                 yggfoi, pch = 19, cex = 2)
## text(p[foi(yggfoi), ], labels = 1:n, col = "white")


## p <- plot2D(yeastgg2[, yeastgg2$condition == "Gly"])
## highlightOnPlot(yeastgg2[, yeastgg2$condition == "Gly"],
##                 yggfoi, pch = 19, cex = 2)
## text(p[foi(yggfoi), ], labels = 1:n, col = "white")


## yyg1 <- yeastgg2[, yeastgg2$replicate == 1]
## yyg2 <- yeastgg2[, yeastgg2$replicate == 2]
## d1 <- pdist(yyg1, "condition")
## d2 <- pdist(yyg2, "condition")

## plot((d1+d2)/2, log2(d1/d2))
## abline(h = 0)
