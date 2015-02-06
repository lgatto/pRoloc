anyUnknown <- function(x, fcol = "markers", unknown = "unknown") 
    any(fData(x)[, fcol] == unknown)

isBinary <- function(x){
    if (inherits(x, "MSnSet"))
        x <- exprs(x)
    x <- as.numeric(x)
    all(unique(x) %in% 0:1)
}

checkFeatureNames <- function(x, y) 
    identical(featureNames(x), featureNames(y))

checkSortedFeatureNames <- function(x, y) 
    identical(sort(featureNames(x)), sort(featureNames(y)))

##' Checks the marker and unknown feature overlap of two \code{MSnSet}
##' instances. 
##'
##' @title Check feature names overlap
##' @param x An \code{MSnSet} instance.
##' @param y An \code{MSnSet} instance.
##' @param fcolx The feature variable to separate unknown
##' (\code{fData(y)$coly == "unknown"}) from the marker features in
##' the \code{x} object.
##' @param fcoly As \code{fcolx}, for the \code{y} object. If missing,
##' the value of \code{fcolx} is used.
##' @param verbose If \code{TRUE} (default), the overlap is printed
##' out on the console.
##' @return Invisibly returns a named list of common markers, unique
##' \code{x} markers, unique \code{y} markers in, common unknowns,
##' unique \code{x} unknowns and unique \code{y} unknowns.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(andy2011)
##' data(andy2011goCC)
##' checkFeatureNamesOverlap(andy2011, andy2011goCC)
##' featureNames(andy2011goCC)[1] <- "ABC"
##' res <- checkFeatureNamesOverlap(andy2011, andy2011goCC)
##' res$markersX
##' res$markersY
checkFeatureNamesOverlap <-
    function(x, y, fcolx = "markers", fcoly,
             verbose = TRUE) {
        if (missing(fcoly)) fcoly <- fcolx
        mx <- featureNames(markerMSnSet(x, fcolx))
        my <- featureNames(markerMSnSet(y, fcoly))
        ux <- featureNames(unknownMSnSet(x, fcoly))
        uy <- featureNames(unknownMSnSet(y, fcoly))
        xym <- intersect(mx, my)
        xm <- setdiff(mx, my)
        ym <- setdiff(my, mx)
        xyu <- intersect(ux, uy)
        xu <- setdiff(ux, uy)
        yu <- setdiff(uy, ux)
        if (verbose) {
            cat("Common markers: ", length(xym), "\n")
            cat("Unique x markers: ", length(xm), "\n")
            cat("Unique y markers: ", length(ym), "\n")
            cat("Common unkowns: ", length(xyu), "\n")
            cat("Unique x unknowns: ", length(xu), "\n")
            cat("Unique y unknowns: ", length(yu), "\n")
        }
        invisible(list(markersXY = xym,
                       markersX = xm,
                       markersY = ym,
                       unknownsXY = xyu,
                       unknownsX = xu,
                       unknownsU = yu))
}
##' Extracts qualitative feature variables from two \code{MSnSet}
##' instances and compares with a contingency table.
##'
##' @title Compare a feature variable overlap
##' @param x An \code{MSnSet} instance. 
##' @param y An \code{MSnSet} instance. 
##' @param fcolx The feature variable to separate unknown
##' (\code{fData(y)$coly == "unknown"}) from the marker features in
##' the \code{x} object.
##' @param fcoly As \code{fcolx}, for the \code{y} object. If missing,
##' the value of \code{fcolx} is used.
##' @param verbose If \code{TRUE} (default), the contingency table of
##' the the feature variables is printed out.
##' @return Invisibly returns a named list with the values of the
##' diagonal, upper and lower triangles of the contingency table.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' res <- checkFvarOverlap(dunkley2006, dunkley2006,
##'                         "markers", "markers.orig")
##' str(res)
checkFvarOverlap <- function(x, y, fcolx = "markers", fcoly,
                                verbose = TRUE) {
    stopifnot(checkFeatureNames(x, y))
    if (missing(fcoly)) fcoly <- fcolx
    mx <- getMarkers(x, fcolx, verbose = FALSE)
    my <- getMarkers(y, fcoly, verbose = FALSE)
    tab <- table(mx, my)
    if (verbose) print(tab)
    invisible(list(matches = diag(tab),
                   lower.mismatches = tab[lower.tri(tab)],
                   upper.mismatches = tab[lower.tri(tab)]))
}
