##' Functions producing a new vector (matrix) marker vector set from
##' an existing matrix (vector) marker set.
##'
##'
##' Sub-cellular markers can be encoded in two different ways. Sets of
##' spatial markers can be represented as character \emph{vectors}
##' (\code{character} or \code{factor}, to be accurate), stored as
##' feature metadata, and proteins of unknown or uncertain
##' localisation (unlabelled, to be classified) are marked with the
##' \code{"unknown"} character. While very handy, this encoding
##' suffers from some drawbacks, in particular the difficulty to label
##' proteins that reside in multiple (possible or actual)
##' localisations. The markers vector feature data is typically named
##' \code{markers}. A new \emph{matrix} encoding is also
##' supported. Each spatial compartment is defined in a column in a
##' binary markers matrix and the resident proteins are encoded with
##' 1s. The markers matrix feature data is typically named
##' \code{Markers}. If proteins are assigned unique localisations only
##' (i.e. no multi-localisation) or their localisation is unknown
##' (unlabelled), then both encodings are equivalent. When the markers
##' are encoded as vectors, features of unknown localisation are
##' defined as \code{fData(object)[, fcol] == "unknown"}. For
##' matrix-encoded markers, unlabelled proteins are defined as
##' \code{rowSums(fData(object)[, fcol]) == 0}.
##'
##' The \code{mrkMatToVec} and \code{mrkVecToMat} functions enable the
##' conversion from matrix (vector) to vector (matrix). The
##' \code{mrkMatAndVec} function generates the missing encoding from
##' the existing one. If the destination encoding already exists, or,
##' more accurately, if the feature variable of the destination
##' encoding exists, an error is thrown. During the conversion from
##' matrix to vector, if multiple possible label exists, they are
##' dropped, i.e. they are converted to \code{"unknown"}. Function
##' \code{isMrkVec} and \code{isMrkMat} can be used to test if a
##' marker set is encoded as a vector or a matrix. \code{mrkEncoding}
##' returns either \code{"vector"} or \code{"matrix"} depending on the
##' nature of the markers.
##' 
##' @title Create a marker vector or matrix.
##' @param object An \code{MSnSet} object 
##' @param vfcol The name of the \emph{vector} marker feature
##' variable. Default is \code{"markers"}.
##' @param mfcol The name of the \emph{matrix} marker feature
##' variable. Default is \code{"Markers"}.
##' @return An updated \code{MSnSet} with a new vector (matrix) marker
##' set.
##' @rdname markers
##' @aliases mrkVecToMat markers
##' @seealso Other functions that operate on markers are
##' \code{\link{getMarkers}}, \code{\link{getMarkerClasses}} and
##' \code{\link{markerMSnSet}}. To add markers to an existing
##' \code{MSnSet}, see the \code{\link{addMarkers}} function and
##' \code{\link{pRolocmarkers}}, for a list of suggested markers.
##' @author Laurent Gatto and Lisa Breckels
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' dunk <- mrkVecToMat(dunkley2006)
##' head(fData(dunk)$Markers)
##' fData(dunk)$markers <- NULL
##' dunk <- mrkMatToVec(dunk)
##' stopifnot(all.equal(fData(dunkley2006)$markers,
##'                     fData(dunk)$markers))
mrkVecToMat <- function(object,
                         vfcol = "markers",
                         mfcol = "Markers") {
    fvl <- fvarLabels(object)
    if (!vfcol %in% fvl) stop(vfcol, " does not exist.")
    if (mfcol %in% fvl) stop(mfcol, " already present.")
    m <- fData(object)[, vfcol]
    um <- levels(factor(m))
    if ("unknown" %in% um) um <- um[um != "unknown"]
    M <- matrix(0, nrow = nrow(object), ncol = length(um))
    rownames(M) <- featureNames(object)
    colnames(M) <- um
    for (j in um) M[which(j == m), j] <- 1
    fData(object)[, mfcol] <- M
    return(object)
}

##' @rdname markers
mrkMatToVec <- function(object,
                        mfcol = "Markers",
                        vfcol = "markers") {
    fvl <- fvarLabels(object)
    if (vfcol %in% fvl) stop(vfcol, " already present.")
    if (!mfcol %in% fvl) stop(mfcol, " does not exist.")
    m <- rep("unknown", nrow(object))
    M <- fData(object)[, mfcol]
    for (i in 1:nrow(object)) {
        k <- M[i, ] == 1
        if (sum(k) == 1) m[i] <- colnames(M)[k]
    }
    fData(object)[, vfcol] <- m
    return(object)
}

##' @rdname markers
mrkMatAndVec <- function(object,
                         vfcol = "markers",
                         mfcol = "Markers") {
    fvl <- fvarLabels(object)
    vex <- vfcol %in% fvl
    mex <- mfcol %in% fvl
    if (!vex & !mex) stop("No marker found.")
    if (vex & !mex)
        object <- mrkVecToMat(object)
    if (!vex & mex)
        object <- mrkMatToVec(object)
    return(object)
}

##' @rdname markers
showMrkMat <- function(object, mfcol = "Markers") {
    M <- fData(object)[, mfcol]
    cat("Localisation count:")
    print(table(rs <- rowSums(M)))
    cat("Single localisations:\n")
    print(apply(M, 2, sum))
    cat("Multiple localisations:")
    multi <- rs > 1
    if (any(multi)) {
        mtab <- table(apply(M[multi, , drop=FALSE],
                    1,
                    function(m)
                        paste0(colnames(M)[m == 1], collapse = "/")))
        print(mtab)
    } else cat("\n  none\n")
}


##' @rdname markers
isMrkMat <- function(object, fcol = "Markers") {
    stopifnot(fcol %in% fvarLabels(object))
    inherits(fData(object)[, fcol], "matrix")
}

##' @rdname markers
##' @param fcol A marker feature variable name.
isMrkVec <- function(object, fcol = "markers") {
    stopifnot(fcol %in% fvarLabels(object))
    isvec <- is.vector(fData(object)[, fcol])
    isfac <- is.factor(fData(object)[, fcol])
    any(c(isvec, isfac))
}

##' @rdname markers
mrkEncoding <- function(object, fcol = "markers") {
    if (isMrkVec(object, fcol)) return("vector")
    if (isMrkMat(object, fcol)) return("matrix")
    stop("Your markers are neither vector nor matrix. See ?markers for details.")
}
