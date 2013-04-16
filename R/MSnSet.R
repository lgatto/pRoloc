##' Convenience accessor to the organelle markers in an 'MSnSet'.
##' This function returns the organelle markers of an
##' \code{MSnSet} instance. As a side effect, it print out a marker table.
##' 
##' @title Returns the organelle markers in an 'MSnSet'
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
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
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The name of the prediction column in the
##' \code{featureData} slot. 
##' @param scol The name of the prediction score column in the
##' \code{featureData} slot. If missing, created by pasting
##' '.scores' after \code{fcol}. If \code{NULL}, ignored.
##' @param t The score threshold. Predictions with score < t are
##' set to 'unknown'. Default is 0.
##' @param verbose If \code{TRUE}, a prediction table is printed and the
##' predictions are returned invisibly. If \code{FALSE}, the predictions
##' are returned.
##' @return A \code{character} of length \code{ncol(object)}. 
##' @author Laurent Gatto
getPredictions <- function(object,
                           fcol,
                           scol,
                           t = 0,
                           verbose = TRUE) {
  stopifnot(!missing(fcol))  
  predictions <- as.character(fData(object)[, fcol])
  if (missing(scol))
    scol <- paste0(fcol, ".scores")
  if (!is.null(scol)) {
    scrs <- fData(object)[, scol]
    predictions[scrs < t] <- "unknown"
  }  
  if (verbose) {
    print(table(predictions))
    invisible(predictions)
  } else {
    return(predictions)
  }
}

##' This functions updates the classification results in an \code{"\linkS4class{MSnSet}"}
##' based on a prediction score threshold \code{t}. All features with a score < t are set
##' to 'unknown'. Note that the original levels are preserved while 'unknown' is added.
##'
##' @title Updates classes based on prediction scores
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The name of the markers column in the \code{featureData} slot. 
##' @param scol The name of the prediction score column in the
##' \code{featureData} slot. If missing, created by pasting
##' '.scores' after \code{fcol}.
##' @param t The score threshold. Predictions with score < t are
##' set to 'unknown'. Default is 0.
##' @return The original \code{object} with a modified \code{fData(object)[, fcol]}
##' feature variable.
##' @author Laurent Gatto
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' ## random scores
##' fData(dunkley2006)$assigned.scores <- runif(nrow(dunkley2006))
##' getPredictions(dunkley2006, fcol = "assigned")
##' getPredictions(dunkley2006, fcol = "assigned", t = 0.5) 
##' x <- minClassScore(dunkley2006, fcol = "assigned", t = 0.5)
##' getPredictions(x, fcol = "assigned")
##' all.equal(getPredictions(dunkley2006, fcol = "assigned", t = 0.5),
##'           getPredictions(x, fcol = "assigned"))
minClassScore <- function(object,
                          fcol,
                          scol,
                          t = 0) {
  stopifnot(!missing(fcol))
  lv <- c(levels(fData(object)[, fcol]),
          "unknown")
  if (missing(scol)) {
    preds <- getPredictions(object, fcol,
                            t = t, verbose = FALSE)
  } else {
    preds <- getPredictions(object, fcol, scol,
                            t = t, verbose = FALSE)
  }
  fData(object)[, fcol] <- factor(preds, levels = lv)
  if (validObject(object))
    object
}

##' This function updates an \code{MSnSet} instances and sets
##' markers class to \code{unknown} if there are less than \code{n}
##' instances. 
##'
##' @title Creates a reduced marker variable
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param n Minumum of marker instances per class.
##' @param fcol The name of the markers column in the \code{featureData}
##' slot. Default is \code{markers}.
##' @return An instance of class \code{"\linkS4class{MSnSet}"} with a new
##' feature variables, named after the original \code{fcol} variable and
##' the \code{n} value. 
##' @author Laurent Gatto
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' d2 <- minMarkers(dunkley2006, 20)
##' getMarkers(dunkley2006)
##' getMarkers(d2, fcol = "markers20")
minMarkers <- function(object, n = 10, fcol = "markers") {
  m <- as.character(fData(object)[, fcol])
  tm <- table(m)
  xx <- names(tm)[tm < n]
  m[m %in% xx] <- "unknown"
  fcol2 <- paste0(fcol, n)
  fData(object)[, fcol2] <- factor(m)
  if (validObject(object))
    return(object)
}


##' The function adds a 'markers' feature variable. These markers are read
##' from a comma separated values (csv) spreadsheet file. This markers file 
##' is expected to have 2 columns (others are ignored) where the first
##' is the name of the marker features and the second the group label.
##' It is essential to assure that  \code{featureNames(object)} and
##' marker names (first column) match, i.e. the same feature identifiers
##' and case fold are used.
##'
##' @title Adds markers to the data
##' @param object An instance of class \code{MSnSet}.
##' @param markerfile A \code{character} with the name the markers' csv file.
##' @param verbose A \code{logical} indicating if number of markers
##' and marker table should be printed to the console.
##' @return A new instance of class \code{MSnSet} with an additional
##' \code{markers} feature variable.
##' @author Laurent Gatto
addMarkers <- function(object, markerfile, verbose = TRUE) {
  if ("markers" %in% fvarLabels(object))
    stop("Detected a feature 'markers' column.")
  mrk <- read.csv(markerfile, stringsAsFactors = FALSE, row.names = 1)
  cmn <- intersect(rownames(mrk), featureNames(object))
  if (length(cmn) == 0) {    
    msg <- paste0("No markers found. Are you sure that the feature names match?\n",
                  "  Feature names: ",
                  paste0(paste(featureNames(object)[1:3], collapse = ", "), "...\n"),
                  "  Markers names: ",
                  paste0(paste(rownames(mrk)[1:3], collapse = ", "), "...\n"))    
    stop(msg)
  }
  if (verbose)
    message("Markers in data: ", length(cmn), " out of ", nrow(object))
  fData(object)$markers <- "unknown"
  fData(object)[cmn, "markers"] <- mrk[cmn, 1]
  object@processingData@processing <-
    c(object@processingData@processing,
      paste0("Added markers from ", basename(markerfile),". ", date()))  
  if (validObject(object)) {
    if (verbose) getMarkers(object)
    return(object)
  }  
}
