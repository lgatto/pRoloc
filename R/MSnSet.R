getRatios <- function(x,log=FALSE) {
  ## x: a vector of numerics
  ## returns a vector of all xi/xj ratios
  x <- as.numeric(x)
  cmb <- combn(length(x),2)
  r <- numeric(ncol(cmb))
  for (i in 1:ncol(cmb)) {
    j <- cmb[1,i]
    k <- cmb[2,i]
    ifelse(log,r[i] <- x[j]-x[k],r[i] <- x[j]/x[k])
  }
  return(r)
}


setGeneric("exprsToRatios",function(object,...) standardGeneric("exprsToRatios"))
setMethod("exprsToRatios",
          "MSnSet",
          function(object,log=FALSE) {            
            r <- apply(exprs(object),1,getRatios,log)
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


##' This function adds an organelle markers column to the \code{featureData}
##' slot of the \code{MSnSet} instance \code{object}. The markers are looked
##' for in an \code{Markers} instance.
##'
##' @title Adds an organelle marker to an \code{MSnSet}
##' @param object An instance of class \code{MSnSet}.
##' @param mrk An instance of class \code{Markers}.
##' @param fcol The name of the \code{featureData} variable (column) to
##' be added. Default is \code{markers}.
##' @return An instance of class \code{MSnSet} with and updated with new 'markers' feature variable.
##' @author Laurent Gatto
addMarkers <- function(object, mrk, fcol = "markers") {
  fData(object)[,fcol] <- "unknown"
  fn <- toupper(featureNames(object)) ## feature names
  mn <- toupper(markers(mrk)$id)      ## marker names
  ## taking visibility of markers into account
  visible <- markers(mrk)$visibility
  mn <- mn[visible]
  cmn <- intersect(mn,fn)
  fn.idx <- match(cmn, fn)
  mn.idx <- match(cmn, mn)
  fData(object)[fn.idx, fcol] <- as.character(markers(mrk)$organelleGOId)[visible][mn.idx]  
  fData(object)[,fcol] <- as.factor(fData(object)[,fcol])
  object@processingData@processing <- c(object@processingData@processing,
                                        paste0("Added markers using ",
                                               MSnbase:::getVariableName(match.call(), "mrk"),
                                               " markers on ", date()))
  if (validObject(object))
    return(object)
}

##' This function removed organelle marker groups that are too small
##' (as defined by \code{n}) to be effectively be used for a complete
##' prediction run, i.e. proper hyperparamter optimisation using
##' stratification and cross-validation and subsequent prediction.
##' The individual protein markers labels are overwritten and set to
##' \code{unknown}.
##'
##' @title Remove small marker groups
##' @param object An instance of class \code{MSnSet}.
##' @param n The minimum accepted group size; marker groups of
##' size <= n will be set to \code{unknown}. Default is 5.
##' @param fcol The feature variable name of the markers. Default
##' is \code{markers}.
##' @param verbose If \code{TRUE} (default), the updated marker 
##' table is printed out.
##' @return An instance of class \code{MSnSet} with updated marker
##' groups.
##' @author Laurent Gatto
minMarkers <- function(object, n = 5,
                       fcol = "markers",
                       verbose = TRUE) {
  ## It would probably be better not to remove these
  ## sparse markers, just disregard them for the class
  ## prediction. Or just set them to unknown to see
  ## if they are wrongly assigned to organelles.
  namesToRemove <- names(which(table(fData(object)[,fcol]) <= n))
  toRemove <- fData(object)[,fcol] %in% namesToRemove
  object <- object[!toRemove,]
  fData(object)[,fcol] <- factor(fData(object)[,fcol])
  if (verbose)
    print(table(fData(object)[,fcol]))
  return(object)
}

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
##' data(dunkley2006)
##' mymarkers <- getMarkers(dunkley2006)
getMarkers <- function(object,
                       fcol = "markers",
                       verbose = TRUE) {
  organelleMarkers <- as.character(fData(object)[, fcol])
  if (verbose) {
    print(table(organelleMarkers))
    return(invisible(organelleMarkers))
  }
  return(organelleMarkers)
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
##' \code{featureData} slot. Default is \code{knn.scores}.
##' @param t The score threshold. Predictions with score < t are
##' set to 'unknown'. Default is 0.
##' @param verbose If \code{TRUE}, a prediction table is printed and the
##' predictions are returned invisibly. If \code{FALSE}, the predictions
##' are returned.
##' @return A \code{character} of length \code{ncol(object)}.
##' @author Laurent Gatto
getPredictions <- function(object,
                           fcol = "knn",
                           scol = "knn.scores",
                           t = 0,
                           verbose = TRUE) {
  predictions <- as.character(fData(object)[, fcol])
  scrs <- fData(object)[, scol]
  predictions[scrs < t] <- "unknown"
  if (verbose)
    print(table(predictions))
  invisible(predictions)
}
