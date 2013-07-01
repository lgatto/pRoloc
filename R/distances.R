setMethod("nndist", "matrix", function(object, ...) nndist_matrix(object))
setMethod("nndist", "MSnSet", function(object, ...) nndist_msnset(object))

.mahalanobis <- function(X, y, S = cov(X)) mahalanobis(X, y, S)

nndist_matrix <- function(object, k = 3,
                   dist = c("euclidian", "mahalanobis"),
                   ...) {
    dist <- match.arg(dist) ## supported distances: euclidian, mahalanobis
    m <- seq(1, k * 2, 2)
    ans <- matrix(NA, nrow = nrow(object), ncol = k * 2)
    rownames(ans) <- rownames(object)
    colnames(ans) <- paste0(rep(c("index", "dist"), k),
                            rep(1:k, each = 2))    
    if (dist == "euclidian") {
        res <- get.knn(object, k = k, ...) # n by k matrix
        ans[, m] <- res$nn.index
        ans[, m + 1] <- res$nn.dist
    } else { ## mahalanobis   
        S <- cov(object)
        res <- apply(object, 1, ## n by n distance matrix
                     function(x) .mahalanobis(object, x, S))
        diag(res) <- Inf 
        idx <- apply(res, 2, order)
        ans[, m] <- t(idx[1:k, ])
        ans[, m + 1] <- matrix(res[idx[1:k, ]], nrow = nrow(object), ncol = k)                
    }
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


