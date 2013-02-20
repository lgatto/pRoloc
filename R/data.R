##' This function takes an instance of class
##' \code{"\linkS4class{MSnSet}"} and sets randomly
##' selected values to NA.
##' The \code{whichNA} can be used to extract the indices
##' of the missing values, as illustrated in the example.
##'
##' @title Create a data with missing values
##' @param object An instance of class \code{MSnSet}.
##' @param nNA The absolute number of missing values to be assigned.
##' @param pNA The proportion of missing values to be assignmed.
##' @return An instance of class \code{MSnSet}, as \code{x},
##' but with the appropriate number/proportion of missing values.
##' The returned object has an additional feature meta-data columns,
##' \code{nNA}
##' @author Laurent Gatto
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' sum(is.na(dunkley2006))
##' dunkleyNA <- makeNaData(dunkley2006, nNA = 150)
##' sum(is.na(dunkleyNA))
##' naIdx <- whichNA(dunkleyNA)
##' head(naIdx)
makeNaData <- function(object, nNA, pNA) {
  stopifnot(inherits(object, "MSnSet"))
  if (missing(nNA) & missing(pNA))
    stop("Provide one of 'nNA' or 'pNA'.")
  if (!missing(nNA) & !missing(pNA))
    stop("Need only one of 'nNA' or 'pNA'.")
  N <- prod(dim(object))
  if (missing(nNA)) {
    if (pNA <= 0 | pNA >= 1)
      stop("Require 0 < pNA < 1")
    nNA <- ceiling(N * pNA)
  }
  if (nNA <= 0 | nNA >= N)
    stop("Require 0 < nNA > ", N)
  .nNA <- sample(N, nNA)
  .naInd <- arrayInd(.nNA, .dim = dim(object))
  .naTab <- table(.naInd[, 1])
  fData(object)$nNA <- 0
  fData(object)$nNA[as.numeric(names(.naTab))] <- .naTab
  exprs(object)[.naInd] <- NA
  stopifnot(sum(fData(object)$nNA) == sum(is.na(object)))
  msg <- {
    if (missing(pNA)) paste0("Set ", nNA, " values to NA")
    else paste0("Set ", nNA, " (", pNA, "%) values to NA")
  }
  msg <- paste(msg, date())  
  object@processingData@processing <-
    c(object@processingData@processing,
      msg)
      
  if (validObject(object))
    return(object)
}

##' @param x A \code{matrix} or an instance of class \code{MSnSet}.
##' @rdname makeNaData
whichNA <- function(x) {
  if (inherits(x, "MSnSet"))
    x <- exprs(x)
  arrayInd(which(is.na(x)), .dim = dim(x))      
}
