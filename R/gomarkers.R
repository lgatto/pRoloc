##' Adds a GO markers to the feature data
##' 
##' @title Add GO markers
##' @param object An instance of class \code{MSnSet}.
##' @param params An instance of class \code{AnnotationParams}. 
##' If missing, \code{\link{getAnnotationParams}} will be used.
##' @param evidence GO evidence filtering.
##' @param useID Logical. Should GO term names or identifiers be used?
##' If \code{TRUE}, identifiers will be used. If \code{FALSE} GO term 
##' names will be used.
##' @param fcol Character. Name of the matrix of markers to be added to the 
##' \code{fData} default is \code{GOMarkers}
##' @param ... Other arguments passed to \code{makeGoSet}
##' @return An updated \code{MSnSet} with new feature data column
##' called \code{GOMarkers} containing a matrix of GO
##' markers
##' @author Lisa M Breckels
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' setAnnotationParams(inputs =
##'                       c("Arabidopsis thaliana genes",
##'                                           "TAIR locus ID"))
##' xx <- addGoMarkers(dunkley2006)
##' dim(fData(xx)$GOMarkers)
##' 
##' ## filter sets
##' xx <- filterMinMarkers(xx, n = 50)
##' dim(fData(xx)$GOMarkers)
##' xx <- filterMaxMarkers(xx, p = .25)
##' dim(fData(xx)$GOMarkers)
##' 
##' ## Subset for specific protein sets
##' sub <- subsetMarkers(xx, keep = c("vacuole"))
##' 
##' ## Order protein sets
##' res <- orderGoMarkers(xx, k = 1:3, p = 1/3, verbose = FALSE)
##' if (interactive()) {
##' pRolocVis(res, fcol = "GOMarkers")
##' }
addGoMarkers <- function(object, params, evidence, 
                         useID = FALSE, fcol = "GOMarkers",
                         ...) {
  if (missing(evidence))
    evidence = NULL
  if (!inherits(params, "AnnotationParams"))
    stop("params must be of class AnnotationParams")
  if (any(fvarLabels(object) == "GOMarkers"))
    stop("colname GOMarkers already exists if fData")
  goSet <- makeGoSet(object, params, evidence = evidence, ...)
  goSet <- filterZeroCols(goSet)
  fData(object)[, fcol] <- exprs(goSet)
  id <- colnames(fData(object)[, fcol])
  orgnames <- goIdToTerm(id, names = FALSE)
  colnames(fData(object)[, fcol]) <- orgnames
  #object <- filterGoMarkers(object) ## Remove any obselete terms
  return (object)
}


##' Removes annotation information that contain less that a
##' certain number/percentage of proteins 
##'
##' @title Removes class/annotation information from a matrix of candidate
##' markers that appear in the \code{fData}.
##' @param object An instance of class \code{MSnSet}.
##' @param n Minimum number of proteins allowed per column. Default is 10.
##' @param p Minimum percentage of proteins per column.
##' @param fcol The name of the matrix of marker information. Default is
##' \code{GOMarkers}.
##' @param verbose Number of marker candidates retained after filtering.
##' @return An updated \code{MSnSet}.
##' @author Lisa M Breckels
##' @seealso \code{addGoMarkers} and example therein.
filterMinMarkers <- function(object, 
                             n = 10, 
                             p,
                             fcol = "GOMarkers",
                             verbose = TRUE) {
  if (!fcol %in% fvarLabels(object)) 
    stop(paste("fcol = ", fcol, "not found in fvarLabels"))
  pm <- fData(object)[, fcol]
  if (is.null(pm)) 
    stop(fcol, " not found.")
  if (!is.matrix(pm)) 
    stop(fcol, " is not a matrix")
  if (!missing(p)) {
    sel <- colSums(pm) > nrow(pm)*p
  } else {
    sel <- colSums(pm) > n
  }
  if (verbose)  
    message("Retaining ", sum(sel), " out of ", ncol(pm), " in ", fcol)
  pm <- pm[, sel, drop = FALSE]
  fData(object)[, fcol] <- pm
  if (validObject(object))
    return(object)
}


##' Removes annotation information that contain more that a
##' certain number/percentage of proteins 
##'
##' @title Removes class/annotation information from a matrix of candidate
##' markers that appear in the \code{fData}.
##' @param object An instance of class \code{MSnSet}.
##' @param n Maximum number of proteins allowed per class/information
##' term.
##' @param p Maximum percentage of proteins per column. Default is 0.2 
##' i.e. remove columns that have information for greater than 20% 
##' of the total number of proteins in the dataset (note: this is useful
##' for example, if information is GO terms, for removing very general
##' and uninformative terms).
##' @param fcol The name of the matrix of marker information. Default is
##' \code{GOMarkers}.
##' @param verbose Number of marker candidates retained after filtering.
##' @return An updated \code{MSnSet}
##' @seealso \code{addGoMarkers} and example therein.
filterMaxMarkers <- function(object, 
                             n, 
                             p = .2,
                             fcol = "GOMarkers",
                             verbose = TRUE) {
  if (!fcol %in% fvarLabels(object)) 
    stop(paste("fcol = ", fcol, "not found in fvarLabels"))
  pm <- fData(object)[, fcol]
  if (is.null(pm)) 
    stop(fcol, " not found.")
  if (!is.matrix(pm)) 
    stop(fcol, " is not a matrix")
  if (!missing(n)) {
    sel <- colSums(pm) < n 
  } else {
    sel <- colSums(pm) < nrow(pm)*p
  }
  if (verbose)  
    message("Retaining ", sum(sel), " out of ", ncol(pm), " in ", fcol)
  pm <- pm[, sel, drop = FALSE]
  fData(object)[, fcol] <- pm
  if (validObject(object))
    return(object)
}



##' Subsets a matrix of markers by specific terms 
##' 
##' @title Subsets markers
##' @param object An instance of class \code{MSnSet}.
##' @param fcol The name of the markers matrix. Default is
##' \code{GOMarkers}.
##' @param keep Integer or character vector specifying the columns to keep 
##' in the markers matrix, as defined by \code{fcol}. 
##' @return An updated \code{MSnSet}
##' @author Lisa M Breckels
##' @seealso \code{addGoMarkers} and example therein.
subsetMarkers <- function(object,
                          fcol = "GOMarkers",
                          keep) {
  if (!fcol %in% fvarLabels(object))
    stop(paste("fcol = ", fcol, "not found in fvarLabels"))
  if (missing(keep))
    stop("No columns selected to subset")
  if (inherits(keep, "integer")) {
    k <- !keep %in% 1:ncol(fData(object)[, fcol])
  } else {
    if (inherits(keep, "character"))
      k <- !keep %in% colnames(fData(object)[, fcol])
    else
      stop("keep must be either integer or character")
  }
  if (sum(k) > 0) {
    missingkeep <- keep[k]
    warning("GO markers ",
            paste(missingkeep, collapse = ", "), " not found")
    keep <- keep[-which(k)]
  }
  .mrkers <- fData(object)[, fcol]
  if (!inherits(.mrkers, "matrix"))
    stop("fcol is not a matrix")
  fData(object)[, fcol] <- .mrkers[, keep, drop = FALSE]
  return(object)
}

##' For a given matrix of candidate markers/annotation information, 
##' this function returns the information ordered according to 
##' the best fit with the data.
##' 
##' As there are typically many protein sets that may fit the data
##' we order protein sets by best fit i.e. cluster tightness, by
##' computing the mean normalised Euclidean distance for all instances 
##' per protein set. 
##' 
##' For each protein set i.e. proteins that have been labelled
##' with a specified term/information criteria, we find the best 
##' \code{k} cluster components for the set (the default is to 
##' test\code{k = 1:5}) according to the minimum mean normalised 
##' pairwise Euclidean distance over all component clusters. 
##' (Note: when testing \code{k} if any components are found to 
##' have less than \code{n} proteins these components are not
##' included and \code{k} is reduced by 1). 
##' 
##' Each component cluster is normalised by \code{N^p} (where 
##' \code{N} is the total number of proteins per component, 
##' and \code{p} is the power). Hueristally, \code{p = 1/3} 
##' and normalising by \code{N^1/3} has been found the optimum 
##' normalisation factor. 
##' 
##' Candidates in the matrix are ordered according to lowest 
##' mean normalised pairwise Euclidean distance as we expect 
##' high density, tight clusters to have the smallest mean 
##' normalised distance.
##' 
##' This function is a wrapper for running \code{clustDist},
##' \code{getNormDist}, see the "Annotating spatial proteomics data"
##' vignette for more details.
##' 
##' @title Orders candidate markers
##' @param object An instance of class \code{MSnSet}.
##' @param fcol The name of the markers matrix. Default is
##' \code{GOMarkers}.
##' @param k The number of clusters to test. Default is \code{k = 1:5} 
##' @param n The minimum number of proteins per component cluster.
##' @param p The normalisation factor, per \code{k} tested
##' @param verbose A \code{logical} indicating if a progress bar should
##' be displayed. Default is \code{TRUE}.
##' @param seed An optional random number generation seed. 
##' @return An updated \code{MSnSet} containing the newly ordered 
##' \code{fcol} matrix.
##' @author Lisa M Breckels
##' @seealso \code{addGoMarkers} and example therein.
orderGoMarkers <- function(object,
                           fcol = "GOMarkers",
                           k = 1:5,
                           n = 5,
                           p = 1/3,
                           verbose = TRUE,
                           seed) {
  ## set seed for reproducibility for kmeans
  if (missing(seed)) {
    seed <- sample(.Machine$integer.max, 1)
  }
  .seed <- as.integer(seed)
  set.seed(.seed)
  ## calculate distances
  message("Calculating GO cluster densities")
  dd <- clustDist(object, k = k, fcol = fcol, n = n, 
                  verbose = verbose, seed = seed)
  ## Normalise by n^1/3
  minDist <- getNormDist(dd, p = 1/3)
  ## Get new order according to lowest distance
  o <- order(minDist)
  ## Re-order `Markers` matrix in `fData`
  fData(object)[, fcol] <- fData(object)[, fcol][, o]
  if (validObject(object)) 
    return(object)
}
