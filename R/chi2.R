setMethod("chi2",
          signature=c(x="numeric",y="numeric"),
          function(x, ## marker intensities - numeric
                   y, ## correlater intensities - numeric
                   method = c("Andersen2003","Wiese2007"),
                   na.rm = FALSE) {
            stopifnot(length(x) == length(y))
            method <- match.arg(method)
            if (length(x)!=length(y))
              stop("marker and correlater must be of same length")
            if (method == "Andersen2003") { ## default when no method specified
              ret <- sum((y-x)^2/length(x), na.rm = na.rm)
            } else {
              ret <- sum((y-x)^2/x, na.rm = na.rm)
            }
            return(ret)
          })

setMethod("chi2",
          signature=c(x="numeric",y="matrix"),
          function(x, ## marker intensities - numeric
                   y, ## correlaters intensities - matrix
                   method = c("Andersen2003","Wiese2007"),
                   na.rm = FALSE) {
            method <- match.arg(method)
            ret <- apply(y,1,function(k) chi2(x, k, method, na.rm))
            return(ret)
          })

setMethod("chi2",
          signature=c(x="matrix",y="numeric"),
          function(x, 
                   y, 
                   method = c("Andersen2003","Wiese2007"),
                   na.rm = FALSE) {
            method <- match.arg(method)
            chi2(y,x,method,na.rm)
          })


setMethod("chi2",
          signature=c(x="matrix",y="matrix"),
          function(x,  ## markers intensities for same organelle - matrix
                   y,  ## correlaters intensities - matrix
                   method = c("Andersen2003","Wiese2007"),
                   fun = NULL, ## function on how to combine chi2 values
                   na.rm = FALSE
                   ) {
            method <- match.arg(method)
            ret <- apply(x, 1, chi2, y, method, na.rm)
            if (is.function(fun))
              ret <- apply(ret,1,fun)
            return(ret)
          })


empPvalues <- function(marker, corMatrix, n=100, ...) {
  ## for 1 marker
  dotargs <- list(...)
  if (is.null(dotargs$na.rm))
    na.rm <- FALSE
  else
    na.rm <- dotargs$na.rm
  method <- dotargs$method
  empP <- rep(NA,nrow(corMatrix))
  chi2res <- chi2(marker, corMatrix, ...)
  for (i in 1:nrow(corMatrix)) {
    rndChi2 <- replicate(n,
                         chi2(sample(marker),
                              corMatrix[i,] ,
                              method, na.rm))
    empP[i] <- sum(chi2res[i] >= rndChi2)
  }
  return(empP/n)
}

