
remap <- function(object, ref = 1) {
    stopifnot(inherits(object, "MSnSetList"))
    x <- msnsets(object)
    refset <- x[[ref]]
    pca1 <- prcomp(exprs(refset), scale = TRUE, center = TRUE)
    colnames(pca1$x) <- sampleNames(refset)
    preds <- lapply(x, function(xx) {
        ans <- predict(pca1, newdata = exprs(xx))
        colnames(ans) <- sampleNames(xx)
        ans
    })
    for (i in seq_along(x)) {
        exprs(object@x[[i]]) <- preds[[i]]
        object@x[[i]] <-
            MSnbase:::logging(object@x[[i]], "remapped (PCA)")
    }
    exprs(object@x[[ref]]) <- pca1$x
    object@log[["remap"]] <- list(type = "PCA", prcomp = pca1, ref = ref)
    if (validObject(object))
        return(object)
}
