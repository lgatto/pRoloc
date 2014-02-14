##' Convenience accessor to the organelle markers in an 'MSnSet'.
##' This function returns the organelle markers of an
##' \code{MSnSet} instance. As a side effect, it print out a marker table.
##' 
##' @title Returns the organelle markers in an 'MSnSet'
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The name of the markers column in the \code{featureData}
##' slot. Default is \code{markers}.
##' @param names A \code{logical} indicating if the markers vector should
##' be named.
##' @param verbose If \code{TRUE}, a marker table is printed and the markers
##' are returned invisibly. If \code{FALSE}, the markers are returned.
##' @return A \code{character} of length \code{ncol(object)}.
##' @author Laurent Gatto
##' @seealso \code{\link{testMarkers}} and \code{\link{minMarkers}}
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' mymarkers <- getMarkers(dunkley2006)
getMarkers <- function(object,
                       fcol = "markers",
                       names = TRUE,
                       verbose = TRUE) {
  organelleMarkers <- as.character(fData(object)[, fcol])
  if (names)
      names(organelleMarkers) <- featureNames(object)
  if (verbose) {
    print(table(organelleMarkers))
    invisible(organelleMarkers)
  } else {
    return(organelleMarkers)
  }
}

##' Tests if the marker class sizes are large enough for the parameter
##' optimisation scheme, i.e. the size is greater that \code{xval + n},
##' where the default \code{xval} is 5 and \code{n} is 2. If the test
##' is unsuccessful, a warning is thrown.
##'
##' In case the test indicates that a class contains too few examples,
##' it is advised to either add some or, if not possible, to remove
##' the class altogether (see \code{\link{minMarkers}})
##' as the parameter optimisation is likely to fail or, at least,
##' produce unreliable results for that class.
##' 
##' @title Tests marker class sizes
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param xval The number cross-validation partitions. See the
##' \code{xval} argument in the parameter optimisation function(s).
##' Default is 5.
##' @param n Number of additional examples. 
##' @param fcol The name of the prediction column in the
##' \code{featureData} slot. Default is \code{"markers"}. 
##' @return If successfull, the test invisibly returns \code{NULL}. Else,
##' it invisibly returns the names of the classes that have too few examples.
##' @author Laurent Gatto
##' @seealso \code{\link{getMarkers}} and \code{\link{minMarkers}}
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' getMarkers(dunkley2006)
##' testMarkers(dunkley2006)
##' toosmall <- testMarkers(dunkley2006, xval = 15)
##' toosmall
testMarkers <- function(object, xval = 5, n = 2, fcol = "markers") {
    mrktab <- table(fData(object)[, fcol])
    N <- xval + 2
    k <- mrktab < N
    ans <- NULL
    if (any(k)) {
        ans <- names(mrktab)[k]
        if (length(ans) == 1) {
            warning(paste(ans, collapse = ", "), " has less than ", N, " markers.")
        } else {
            warning(paste(ans, collapse = ", "), " have/has less than ", N, " markers.")
        }
    } 
    invisible(ans)
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
##' and case fold are used. Alternatively, a markers named vector as provided
##' by the \code{\link{pRolocmarkers}} function can also be used.
##'
##' @title Adds markers to the data
##' @param object An instance of class \code{MSnSet}.
##' @param markers A \code{character} with the name the markers' csv
##' file or a named character of markers as provided by \code{\link{pRolocmarkers}}. 
##' @param verbose A \code{logical} indicating if number of markers
##' and marker table should be printed to the console.
##' @return A new instance of class \code{MSnSet} with an additional
##' \code{markers} feature variable.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' atha <- pRolocmarkers("atha")
##' try(addMarkers(dunkley2006, atha)) ## markers already exists
##' fData(dunkley2006)$markers.org <- fData(dunkley2006)$markers
##' fData(dunkley2006)$markers <- NULL
##' marked <- addMarkers(dunkley2006, atha)
##' fvarLabels(marked)
##' plot2D(marked)
##' addLegend(marked, where = "topleft", cex = .7)
addMarkers <- function(object, markers, verbose = TRUE) {
  if ("markers" %in% fvarLabels(object))
    stop("Detected an existing 'markers' feature column.")
  if (length(markers) == 1 && file.exists(markers)) {
      mrk <- read.csv(markers, stringsAsFactors = FALSE, row.names = 1)
      mfrom <- basename(markers)
  } else {
      mrk <- cbind(markers)
      mfrom <- paste0(" '",
                      MSnbase:::getVariableName(match.call(), "markers"),
                      "' marker vector")
  }
  dups <- duplicated(rownames(mrk))
  if (any(dups))
      stop("Please remove duplicated entries in your markers:",
           paste(rownames(mrk)[dups], collapse = " "))
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
      paste0("Added markers from ", mfrom,". ", date()))  
  if (validObject(object)) {
    if (verbose) getMarkers(object)
    return(object)
  }  
}

##' This function extracts the marker proteins into a new
##' \code{MSnSet}. 
##' 
##' @title Extract marker subset
##' @param object An instance of class \code{MSnSet}
##' @param fcol The name of the feature data column, that
##' will be used to separate the markers from the proteins
##' of unknown localisation (with
##' \code{fData(object)[, fcol] == "unknown")}).
##' Default is to use \code{"markers"}.
##' @return An new \code{MSnSet} with marker/unknown proteins only.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' mrk <- markerSet(dunkley2006)
##' unk <- unknownSet(dunkley2006)
##' dim(dunkley2006)
##' dim(mrk)
##' dim(unk)
##' table(fData(dunkley2006)$markers)
##' table(fData(mrk)$markers)
##' table(fData(unk)$markers)
markerSet <- function(object, fcol = "markers") {
  mrk <- fData(object)[, fcol]
  object <- object[mrk != "unknown", ]
  ## drop "unknown" level
  fData(object)[, fcol] <- factor(fData(object)[, fcol])
  if (validObject(object))
    return(object)
}

##' @rdname markerSet
unknownSet <- function(object, fcol = "markers") {
  mrk <- fData(object)[, fcol]
  object <- object[mrk == "unknown", ]
  fData(object)[, fcol] <- factor(fData(object)[, fcol])
  if (validObject(object))
    return(object)
}
