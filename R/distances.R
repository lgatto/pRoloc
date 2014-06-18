setMethod("nndist", "matrix",
          function(object, k = 3, dist = c("euclidean", "mahalanobis"), ...)
          nndist_matrix(object, k, match.arg(dist), ...))
setMethod("nndist", "MSnSet",
          function(object, k = 3, dist = c("euclidean", "mahalanobis"), ...)
          nndist_msnset(object, k, match.arg(dist), ...))

.mahalanobis <- function(X, y, S = cov(X)) mahalanobis(X, y, S)

nndist_matrix <- function(object, k, dist = dist, ...) {
    m <- seq(1, k * 2, 2)
    ans <- matrix(NA, nrow = nrow(object), ncol = k * 2)
    rownames(ans) <- rownames(object)
    colnames(ans) <- paste0(rep(c("index", "dist"), k),
                            rep(1:k, each = 2)) 
    if (dist == "euclidean") {
        res <- get.knn(object, k = k, ...) # n by k matrix
        ans[, m] <- res$nn.index
        ans[, m + 1] <- res$nn.dist
        colnames(ans) <- paste0(colnames(ans), "euc")
    } else { ## mahalanobis   
        S <- cov(object)
        res <- apply(object, 1, ## n by n distance matrix
                     function(x) .mahalanobis(object, x, S))
        diag(res) <- Inf 
        idx <- apply(res, 2, order)
        ans[, m] <- t(idx[1:k, ])
        ans[, m + 1] <- matrix(res[idx[1:k, ]], nrow = nrow(object), ncol = k)
        colnames(ans) <- paste0(colnames(ans), "mah")
    }
    return(ans)
}

nndist_msnset <- function(object, k, dist = dist, ...) {
    res <- nndist_matrix(exprs(object), k, dist, ...)                  
    fData(object) <- cbind(fData(object),
                           res)
    if (validObject(object))
        return(object)
}

## Use of 'get.knnx' instead of 'get.knn' allows one to have a query matrix
nndistx_matrix <- function(object, query, k, ...) {
  m <- seq(1, k * 2, 2)
  ans <- matrix(NA, nrow = nrow(query), ncol = k * 2)
  rownames(ans) <- rownames(query)
  colnames(ans) <- paste0(rep(c("index", "dist"), k),
                          rep(1:k, each = 2)) 
  res <- get.knnx(object, query, k = k, ...) # n by k matrix
  ans[, m] <- res$nn.index
  ans[, m + 1] <- res$nn.dist
  colnames(ans) <- paste0(colnames(ans), "euc")
  return(ans)
}

