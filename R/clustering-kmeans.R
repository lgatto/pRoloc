.optim <- function(k, object, cl,
                   criterion = c("BIC", "AIC"), ...) {
  
  criterion <- match.arg(criterion, several.ok = TRUE) 
  .var <- function(x) 
    sum(x^2 - mean(x)^2)/length(x)
  
  X <- exprs(object)
  N <- nrow(X)
  P <- ncol(X)
  if (missing(cl))
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
  if (length(criterion) == 2) {
    ans <- list(BIC = .BIC, AIC = .AIC)
  } else {
    ans <- switch(criterion,
                  BIC = .BIC,
                  AIC = .AIC)
  }
  return(ans)
}


setGeneric("kmeansClustering",
           function(object, params, ...)
           standardGeneric("kmeansClustering"))

setMethod("kmeansClustering",
          signature(object = "MSnSet",
                    params = "missing"),
          function(object, ...) {
            cl <- kmeans(exprs(object), ...)
            fData(object)$kmeans <- cl$cluster
            ## possibly update processingData
            if (validObject(object))
              return(object)
          })

setMethod("kmeansClustering",
          signature(object = "matrix",
                    params = "missing"),
          function(object, ...) {
            kmeans(object, ...)
          })

setMethod("kmeansClustering",
          signature(object = "MSnSet",
                    params = "ClustRegRes"),
          function(object, params){
            ## extract params
            ## object <- kmeans(object, params)
            ## return(object)          
          })

setMethod("kmeansOptimisation",
          signature(object = "MSnSet"),
          function(object, k = 1:20, ...) {
            sapply(1:20, .optim, object, ...)
          })
