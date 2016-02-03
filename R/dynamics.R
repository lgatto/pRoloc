
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
    if (!MSnbase:::.sameNbCol(object))
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
