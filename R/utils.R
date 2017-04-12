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
                       unknownsY = yu))
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


##' This function replaces a string or regular expression in a feature
##' variable using the \code{\link{sub}} function.
##'
##' @title Update a feature variable
##' @param object An instance of class \code{MSnSet}.
##' @param fcol Feature variable to be modified. Default is
##'     \code{"markers"}. If \code{NULL}, all feature variables will
##'     updated.
##' @param from A \code{character} defining the string or regular
##'     expression of the pattern to be replaced. Default is the empty
##'     string, i.e. the regular expression \code{"^$"}.  See
##'     \code{\link{sub}} for details. If \code{NA}, then \code{NA}
##'     values are replaced by \code{to}.
##' @param to A replacement for matched pattern. Default is
##'     \code{"unknown"}.  See \code{\link{sub}} for details.
##' @param ... Additional arguments passed to \code{\link{sub}}.
##' @return An updated \code{MSnSet}.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' getMarkers(dunkley2006, "markers")
##' dunkley2006 <- fDataToUnknown(dunkley2006,
##'                               from = "unknown", to = "unassigned")
##' getMarkers(dunkley2006, "markers")
fDataToUnknown <- function(object, fcol = "markers",
                           from = "^$", to = "unknown",
                           ...) {
    if (is.null(fcol))
        fcol <- fvarLabels(object)
    for (.fcol in fcol) {
        if (is.na(from))
            fData(object)[, .fcol][is.na(fData(object)[, .fcol])] <- to
        else 
            fData(object)[, .fcol] <-
                sub(from, to, fData(object)[, .fcol], ...)
    }
    if (validObject(object))
        return(object)
}

##' Calculates class weights to be used for parameter optimisation and
##' classification such as \code{\link{svmOptimisation}} or
##' \code{\link{svmClassification}} - see the \emph{pRoloc tutorial}
##' vignette for an example. The weights are calculated for all
##' non-\emph{unknown} classes the inverse of the number of
##' observations.
##'
##' @title Calculate class weights
##' @param object An instance of class \code{MSnSet}
##' @param fcol The name of the features to be weighted
##' @return A \code{table} of class weights
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(hyperLOPIT2015)
##' classWeights(hyperLOPIT2015)
##' data(dunkley2006)
##' classWeights(dunkley2006)
classWeights <- function(object, fcol = "markers") {
    stopifnot(inherits(object, "MSnSet"))
    w <- table(fData(object)[, fcol])
    w <- 1/w[names(w) != "unknown"]
    w
}

##' A function to calculate average marker profiles.
##'
##' @title Marker consensus profiles
##' @param object An instance of class \code{MSnSet}.
##' @param fcol Feature meta-data label (fData column name) defining
##'     the groups to be differentiated using different
##'     colours. Default is \code{markers}.
##' @param method A \code{function} to average marker
##'     profiles. Default is \code{mean}.
##' @return A \code{matrix} of dimensions \emph{number of clusters}
##'     (exluding unknowns) by \emph{number of fractions}.
##' @author Laurent Gatto and Lisa M. Breckels
##' @seealso The \code{\link{mrkHClust}} function to produce a
##'     hierarchical cluster.
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' mrkConsProfiles(dunkley2006)
##' mrkConsProfiles(dunkley2006, method = median)
##' mm <- mrkConsProfiles(dunkley2006)
##' ## Reorder fractions
##' o <- order(dunkley2006$fraction)
##' ## Plot mean organelle profiles using the
##' ## default pRoloc colour palette.
##' matplot(t(mm[, o]), type = "l",
##'         xlab = "Fractions", ylab = "Relative intensity",
##'         main = "Mean organelle profiles",
##'         col = getStockcol(), lwd = 2, lty = 1)
##' ## Add a legend
##' addLegend(markerMSnSet(dunkley2006), where = "topleft")
mrkConsProfiles <- function(object, fcol = "markers", method = mean) {
    object <- markerMSnSet(object, fcol)
    cl <- getMarkerClasses(object, fcol)
    ind <- lapply(cl, function(x) which(fData(object)[, fcol] == x))
    names(ind) <- cl
    profs <- lapply(ind, function(x) exprs(object)[x, , drop = FALSE])
    mm <- t(sapply(profs, function(z) apply(z, 2, method)))
    return(mm)
}
