.optim <- function(k, X, 
                   criterion = c("BIC", "AIC"),
                   ...) {
  ## value must a named vector
  criterion <- match.arg(criterion, several.ok = TRUE) 
  .var <- function(x) 
    sum(x^2 - mean(x)^2)/length(x)
  
  N <- nrow(X)
  P <- ncol(X)
  cl <- kmeans(X, centers = k, ...)
  Y <- cl$cluster
  Nc <- cl$size
  K <- length(Nc)

  ## per cluster variance - using N instead of N-1 for clusters of size 1
  Vc <- matrix(0, nrow = P, ncol = K)
  rownames(Vc) <- colnames(X)
  colnames(Vc) <- 1:K
  for (k in 1:K) {
    for (p in 1:P) {
      Vc[p, k] <- .var(X[Y == k, p]) 
    }
  }
  ## total variance
  V <- matrix(apply(X, 2, var), ncol = K, nrow = P)
  LL <- -Nc * colSums(log(Vc + V)/2) 
  .BIC <- -2 * sum(LL) + 2 * K * P * log(N)
  .AIC <- -2 * sum(LL) + 4 * K * P
  ans <- c(BIC = .BIC, AIC = .AIC)
  if (length(criterion) == 1) 
    ans <- switch(criterion,
                  BIC = ans["BIC"],
                  AIC = ans["AIC"])
  return(ans)
}

## setMethod("kmeansClustering",
##           signature(object = "MSnSet",
##                     params = "missing"),
##           function(object, ...) {
##             cl <- kmeans(exprs(object), ...)
##             fData(object)$kmeans <- cl$cluster
##             ## possibly update processingData
##             if (validObject(object))
##               return(object)
##           })

## setMethod("kmeansClustering",
##           signature(object = "matrix",
##                     params = "missing"),
##           function(object, ...) {
##             kmeans(object, ...)
##           })

## setMethod("kmeansClustering",
##           signature(object = "MSnSet",
##                     params = "ClustRegRes"), 
##           function(object, params, criterion = "BIC"){
##             k <- getParams(params, criterion)
##             if (length(params@algoparams) == 0) {
##               cl <- kmeans(exprs(object),
##                            centers = getParams(params, criterion))
##             } else {
##               cl <- kmeans(exprs(object),
##                            centers = getParams(params, criterion),
##                            as.pairlist(params@algoparams))
##             }
##             fData(object)$kmeans <- cl$cluster
##             if (validObject(object))
##               return(object)
##           })


## setMethod("kmeansOptimisation",
##           signature(object = "MSnSet"), 
##           function(object, k = 1:20, ...) {
##             x <- lapply(k, .optim, exprs(object), ...)
##             .criteria <- names(x[[1]])
##             x <- do.call(rbind, x)
##             x <- cbind(k, x)
##             ans <- new("ClustRegRes",
##                        algorithm = "kmeans",
##                        criteria = .criteria,
##                        parameters = list(k = k),
##                        results = x,
##                        algoparams = list(...))
##             if (validObject(ans))
##               return(ans)
##           })

## - compare/contrast 2 clusters - could be 2 fcols
## - estimate parameters - use 1 fcol

## setGeneric("kmeansOptimisation",
##             function(object, fcol, ...) standardGeneric("kmeansOptimisation"))
## setMethod("kmeansOptimisation",
##           signature(object = "MSnSet",
##                     fcol = "character"),
##           function(object, model, ...) {
##             ## compare results of cl and fData(object)[, fcol]
##             return(NULL)
##           })
