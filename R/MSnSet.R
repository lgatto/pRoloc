getRatios <- function(x, log = FALSE) {
  ## x: a vector of numerics
  ## returns a vector of all xi/xj ratios
  x <- as.numeric(x)
  cmb <- combn(length(x),2)
  r <- numeric(ncol(cmb))
  for (i in 1:ncol(cmb)) {
    j <- cmb[1,i]
    k <- cmb[2,i]
    ifelse(log,
           r[i] <- x[j]-x[k],
           r[i] <- x[j]/x[k])
  }
  return(r)
}


setMethod("exprsToRatios",
          "MSnSet",
          function(object, log = FALSE) {            
            r <- apply(exprs(object), 1, getRatios, log)
            r <- t(r)
            rownames(r) <- featureNames(object)
            cmb <- combn(ncol(object),2)            
            ratio.description <- apply(cmb,2, function(x)
                                       paste(sampleNames(object)[x[1]],
                                             sampleNames(object)[x[2]],
                                             sep="/"))
            phenodata <- new("AnnotatedDataFrame",
                             data=data.frame(ratio.description))
            processingdata <- processingData(object)
            processingdata@processing <- c(processingdata@processing,
                                           paste("Intensities to ratios: ",date(),sep=""))
            message("Dropping protocolData.")
            res <- new("MSnSet",
                       exprs = r,
                       featureData = featureData(object),
                       phenoData = phenodata,
                       processingData = processingdata,
                       experimentData = experimentData(object))
            if (validObject(res))
              return(res)
          })



##' Convenience accessor to the organelle markers in an 'MSnSet'.
##' This function returns the organelle markers of an
##' \code{MSnSet} instance. As a side effect, it print out a marker table.
##' 
##' @title Returns the organelle markers in an 'MSnSet'
##' @param object An instance of class \code{MSnSet}.
##' @param fcol The name of the markers column in the \code{featureData}
##' slot. Default is \code{markers}.
##' @param verbose If \code{TRUE}, a marker table is printed and the markers
##' are returned invisibly. If \code{FALSE}, the markers are returned.
##' @return A \code{character} of length \code{ncol(object)}.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' mymarkers <- getMarkers(dunkley2006)
getMarkers <- function(object,
                       fcol = "markers",
                       verbose = TRUE) {
  organelleMarkers <- as.character(fData(object)[, fcol])
  if (verbose) {
    print(table(organelleMarkers))
    invisible(organelleMarkers)
  } else {
    return(organelleMarkers)
  }
}


##' Convenience accessor to the predicted feature localisation in an 'MSnSet'.
##' This function returns the predictions of an
##' \code{MSnSet} instance. As a side effect, it prints out a prediction table.
##' 
##' @title Returns the predictions in an 'MSnSet'
##' @param object An instance of class \code{MSnSet}.
##' @param fcol The name of the prediction column in the
##' \code{featureData} slot. Default is \code{knn}.
##' @param scol The name of the prediction score column in the
##' \code{featureData} slot. If missing, created by pasting
##' '.scores' after \code{fcol}.
##' @param t The score threshold. Predictions with score < t are
##' set to 'unknown'. Default is 0.
##' @param verbose If \code{TRUE}, a prediction table is printed and the
##' predictions are returned invisibly. If \code{FALSE}, the predictions
##' are returned.
##' @return A \code{character} of length \code{ncol(object)}. 
##' @author Laurent Gatto
getPredictions <- function(object,
                           fcol = "knn",
                           scol,
                           t = 0,
                           verbose = TRUE) {
  predictions <- as.character(fData(object)[, fcol])
  if (missing(scol))
    scol <- paste0(fcol, ".scores")
  scrs <- fData(object)[, scol]
  predictions[scrs < t] <- "unknown"
  if (verbose) {
    print(table(predictions))
    invisible(predictions)
  } else {
    return(predictions)
  }
}
