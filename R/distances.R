setMethod("nndist", "matrix", function(object, ...) nndist_matrix(object))
setMethod("nndist", "MSnSet", function(object, ...) nndist_msnset(object))

nndist_matrix <- function(object, k = 3,
                   dist = "euclidian" ,
                   ...) {
    dist <- match.arg(dist)
    ## supported distances: euclidian
    ## to do: mahalanobis, possibly other from dist
    res <- get.knn(object, k = k, ...)
    m <- seq(1, k * 2, 2)
    ans <- matrix(NA, nrow = nrow(object), ncol = k * 2)
    ans[, m] <- res$nn.index
    ans[, m + 1] <- res$nn.dist
    rownames(ans) <- rownames(object)
    colnames(ans) <- paste0(rep(c("index", "dist"), k),
                            rep(1:k, each = 2))
    return(ans)
}

nndist_msnset <- function(object, k = 3,
                          dist = "euclidian" ,
                          ...) {
    res <- nndist(exprs(object), k = k,
                  dist = dist, ...)
    fData(object) <- cbind(fData(object),
                           res)
    if (validObject(object))
        return(object)
}

    
