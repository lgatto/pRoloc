##' Convenience accessor to the organelle markers in an \code{MSnSet}.
##' This function returns the organelle markers of an \code{MSnSet}
##' instance. As a side effect, it print out a marker table.
##' 
##' @title Get the organelle markers in an \code{MSnSet}
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The name of the markers column in the \code{featureData}
##' slot. Default is \code{"markers"}.
##' @param names A \code{logical} indicating if the markers vector should
##' be named. Ignored if markers are encoded as a matrix.
##' @param verbose If \code{TRUE}, a marker table is printed and the markers
##' are returned invisibly. If \code{FALSE}, the markers are returned.
##' @return A \code{character} (\code{matrix}) of length (ncol)
##' \code{ncol(object)}, depending on the vector or matrix encoding of
##' the markers. 
##' @author Laurent Gatto
##' @seealso See \code{\link{getMarkerClasses}} to get the classes
##' only. See \code{\link{markers}} for details about spatial markers
##' storage and encoding.
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' ## marker vectors
##' myVmarkers <- getMarkers(dunkley2006)
##' head(myVmarkers)
##' ## marker matrix
##' dunkley2006 <- mrkVecToMat(dunkley2006, mfcol = "Markers")
##' myMmarkers <- getMarkers(dunkley2006, fcol = "Markers")
##' head(myMmarkers)
getMarkers <- function(object,
                       fcol = "markers",
                       names = TRUE,
                       verbose = TRUE)
    switch(mrkEncoding(object, fcol),
           vector = getVecMarkers(object, fcol, names, verbose),
           matrix = getMatMarkers(object, fcol, verbose))

getVecMarkers <- function(object, fcol, names, verbose) {
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

getMatMarkers <- function(object, fcol, verbose) {
    if (verbose) {
        showMrkMat(object, fcol)
        invisible(fData(object)[, fcol])
    }
    return(fData(object)[, fcol])
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
##' @param error A \code{logical} specifying if an error should be
##' thown, instead of a warning.
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
##' try(testMarkers(dunkley2006, xval = 15, error = TRUE))
testMarkers <- function(object, xval = 5, n = 2,
                        fcol = "markers", error = FALSE) {
    mrktab <- table(fData(object)[, fcol])
    N <- xval + 2
    k <- mrktab < N
    ans <- NULL
    if (any(k)) {
        ans <- names(mrktab)[k]
        if (length(ans) == 1) {
            msg <- paste0(paste(ans, collapse = ", "),
                          " has less than ", N, " markers.")
        } else {
            msg <- paste0(paste(ans, collapse = ", "),
                          " have/has less than ", N, " markers.")
        }
        if (error) stop(msg)
        else warning(msg)
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
##' '.scores' after \code{fcol}. 
##' @param mcol The feature meta data column containing the labelled training 
##' data.
##' @param t The score threshold. Predictions with score < t are set
##' to 'unknown'. Default is 0. It is also possible to define
##' thresholds for each prediction class, in which case, \code{t} is a
##' named numeric with names exactly matching the unique prediction
##' class names.
##' @param verbose If \code{TRUE}, a prediction table is printed and the
##' predictions are returned invisibly. If \code{FALSE}, the predictions
##' are returned.
##' @return An instance of class "\linkS4class{MSnSet}" with \code{fcol.pred} feature
##' variable storing the prediction results according to the chosen threshold.
##' @author Laurent Gatto and Lisa Breckels
##' @seealso \code{\link{orgQuants}} for calculating organelle-specific
##' thresholds.
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' res <- svmClassification(dunkley2006, fcol = "pd.markers",
##'                          sigma = 0.1, cost = 0.5)
##' fData(res)$svm[500:510]
##' fData(res)$svm.scores[500:510]
##' getPredictions(res, fcol = "svm", t = 0) ## all predictions
##' getPredictions(res, fcol = "svm", t = .9) ## single threshold 
##' ## 50% top predictions per class
##' ts <- orgQuants(res, fcol = "svm", t = .5)
##' getPredictions(res, fcol = "svm", t = ts)
getPredictions <- function(object,
                           fcol,
                           scol,
                           mcol = "markers",
                           t = 0,
                           verbose = TRUE) {
  stopifnot(!missing(fcol))
    if (missing(scol))
        scol <- paste0(fcol, ".scores")
    ans <- predictions <-
        as.character(fData(object)[, fcol])
    predclasses <- unique(predictions)
    ## Note: If any of the thresholds are NA we set them to Infinity so  
    ## no new assignments can be made. Usually, a threshold is calculated
    ## from the distribution of class scores of the unlabelled data.
    ## However, if there are no new assignments for a particular 
    ## class, there are no scores on which to calculate the threshold 
    ## and this can result in NA values.
    if (anyNA(t)) {
      t[whichNA(t)] <- Inf
      warning('t contains NA, setting t to Inf')
    }
    if (length(t) > 1) {
        if (!all(sort(names(t)) == sort(predclasses)))
            stop("Class-specific score names do not match the class names exactly:\n",
                 "   score names: ", paste(sort(names(t)), collapse = ", "), "\n",
                 "   class names: ", paste(sort(predclasses), collapse = ", "))
        tt <- as.vector(t[predictions])
        ans <- ifelse(fData(object)[, scol] < tt,
                      "unknown", predictions)
    } else {
        scrs <- fData(object)[, scol]
        ans[scrs < t] <- "unknown"
    }
    train <- as.character(fData(object)[, mcol])
    train.ind <- which(train != "unknown")
    ans[train.ind] <- train[train.ind]
    t <- format(t, digits = 2)
    if (length(t) > 1)
      p <- paste("thresholds:", paste(paste(names(t), t, sep = " = "), collapse = ", "))
    else
      p <- paste("global threshold =", t)
    l <- paste0(fcol, ".pred")
    fData(object)[, l] <- ans
    if (verbose) {
        print(table(ans))
        invisible(ans)
    }
    object@processingData@processing <- c(processingData(object)@processing, 
                                        paste("Added", fcol, "predictions according to", p, date()))
    return(object)
}


minClassScore <- function(object,
                          fcol,
                          scol,
                          t = 0) {
    .Deprecated("getPredictions")
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
##' @seealso \code{\link{getPredictions}} to filter based on
##' classification scores.
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


##' The function adds a 'markers' feature variable. These markers are
##' read from a comma separated values (csv) spreadsheet file. This
##' markers file is expected to have 2 columns (others are ignored)
##' where the first is the name of the marker features and the second
##' the group label. Alternatively, a markers named vector as provided
##' by the \code{\link{pRolocmarkers}} function can also be used.
##'
##' It is essential to assure that \code{featureNames(object)} (or
##' \code{fcol}, see below) and marker names (first column) match,
##' i.e. the same feature identifiers and case fold are used.
##'
##' @title Adds markers to the data
##' @param object An instance of class \code{MSnSet}.
##' @param markers A \code{character} with the name the markers' csv
##' file or a named character of markers as provided by
##' \code{\link{pRolocmarkers}}.
##' @param mcol A \code{character} of length 1 defining the feature
##' variable label for the newly added markers. Default is
##' \code{"markers"}.
##' @param fcol An optional feature variable to be used to match
##' against the markers. If missing, the feature names are used.
##' @param verbose A \code{logical} indicating if number of markers
##' and marker table should be printed to the console.
##' @return A new instance of class \code{MSnSet} with an additional
##' \code{markers} feature variable.
##' @seealso See \code{\link{pRolocmarkers}} for a list of spatial
##' markers and \code{\link{markers}} for details about markers
##' encoding.
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
##' ## if 'makers' already exists
##' marked <- addMarkers(marked, atha, mcol = "markers2")
##' fvarLabels(marked)
##' stopifnot(all.equal(fData(marked)$markers, fData(marked)$markers2))
##' plot2D(marked)
##' addLegend(marked, where = "topleft", cex = .7)
addMarkers <- function(object, markers,
                       mcol = "markers",
                       fcol, verbose = TRUE) {
    if (mcol %in% fvarLabels(object))
        stop("Detected an existing '", mcol, "' feature column.")
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
    if (missing(fcol)) {
        fn <- featureNames(object)
    } else {
        if (!fcol %in% fvarLabels(object))
            stop("'", fcol, "' not found in feature variables.")
        fn <- as.character(fData(object)[, fcol])
    }
    cmn <- fn %in% rownames(mrk)

    if (sum(cmn) == 0) {
        msg <- paste0("No markers found. Are you sure that the feature names match?\n",
                      "  Feature names: ",
                      paste0(paste(featureNames(object)[1:3], collapse = ", "), "...\n"),
                      "  Markers names: ",
                      paste0(paste(rownames(mrk)[1:3], collapse = ", "), "...\n"))
        stop(msg)
    }
    if (verbose)
        message("Markers in data: ", sum(cmn), " out of ", nrow(object))
    k <- match(fn[cmn], rownames(mrk))
    fData(object)[, mcol] <- "unknown"
    fData(object)[cmn, mcol] <- mrk[k, 1]
    object@processingData@processing <-
        c(object@processingData@processing,
          paste0("Added markers from ", mfrom,". ", date()))
    if (validObject(object)) {
        if (verbose) getMarkers(object, fcol = mcol)
        return(object)
    }
}

##' These function extract the marker or unknown proteins into a new
##' \code{MSnSet}.
##' 
##' @title Extract marker/unknown subsets
##' @param object An instance of class \code{MSnSet}
##' @param fcol The name of the feature data column, that will be used
##' to separate the markers from the proteins of unknown
##' localisation. When the markers are encoded as vectors, features of
##' unknown localisation are defined as \code{fData(object)[, fcol] ==
##' "unknown"}. For matrix-encoded markers, unlabelled proteins are
##' defined as \code{rowSums(fData(object)[, fcol]) == 0}. Default is
##' \code{"markers"}.
##' @return An new \code{MSnSet} with marker/unknown proteins only.
##' @seealso \code{\link{sampleMSnSet}} \code{\link{testMSnSet}} and
##' \code{\link{markers}} for markers encoding.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' mrk <- markerMSnSet(dunkley2006)
##' unk <- unknownMSnSet(dunkley2006)
##' dim(dunkley2006)
##' dim(mrk)
##' dim(unk)
##' table(fData(dunkley2006)$markers)
##' table(fData(mrk)$markers)
##' table(fData(unk)$markers)
##' ## matrix-encoded markers
##' dunkley2006 <- mrkVecToMat(dunkley2006)
##' dim(markerMSnSet(dunkley2006, "Markers"))
##' stopifnot(all.equal(featureNames(markerMSnSet(dunkley2006, "Markers")),
##'                     featureNames(markerMSnSet(dunkley2006, "markers"))))
##' dim(unknownMSnSet(dunkley2006, "Markers"))
##' stopifnot(all.equal(featureNames(unknownMSnSet(dunkley2006, "Markers")),
##'                     featureNames(unknownMSnSet(dunkley2006, "markers"))))
markerMSnSet <- function(object, fcol = "markers")
    switch(mrkEncoding(object, fcol),
           vector = vecMarkerMSnSet(object, fcol),
           matrix = matMarkerMSnSet(object, fcol))

vecMarkerMSnSet <- function(object, fcol) {
    mrk <- fData(object)[, fcol]
    object <- object[mrk != "unknown", ]
    ## drop "unknown" level
    fData(object)[, fcol] <- factor(fData(object)[, fcol])
    if (validObject(object))
        return(object)
}

matMarkerMSnSet <- function(object, fcol) {
    rs <- rowSums(fData(object)[, fcol])
    object <- object[rs > 0, ]
    if (validObject(object))
        return(object)
}

##' @rdname markerMSnSet
unknownMSnSet <- function(object, fcol = "markers")
    switch(mrkEncoding(object, fcol),
           vector = vecUnknownMSnSet(object, fcol),
           matrix = matUnknownMSnSet(object, fcol))

vecUnknownMSnSet <- function(object, fcol) {
    mrk <- fData(object)[, fcol]
    object <- object[mrk == "unknown", ]
    fData(object)[, fcol] <- factor(fData(object)[, fcol])
    if (validObject(object))
        return(object)
}

matUnknownMSnSet <- function(object, fcol) {
    rs <- rowSums(fData(object)[, fcol])
    object <- object[rs == 0, ]
    if (validObject(object))
        return(object)
}


##' This function creates a stratified 'test' \code{MSnSet} which can be used 
##' for algorihtmic development. A \code{"\linkS4class{MSnSet}"} containing only
##' the marker proteins, as defined in \code{fcol}, is returned with a new 
##' feature data column appended called \code{test} in which a stratified subset
##' of these markers has been relabelled as 'unknowns'.
##' 
##' @title Create a stratified 'test' \code{MSnSet}
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}
##' @param fcol The feature meta-data column name containing the
##' marker definitions on which the data will be stratified. Default
##' is \code{markers}.
##' @param size The size of the data set to be extracted. Default is
##' 0.2 (20 percent).
##' @param seed The optional random number generator seed.
##' @return An instance of class \code{"\linkS4class{MSnSet}"} which
##' contains only the proteins that have a labelled localisation
##' i.e. the marker proteins, as defined in \code{fcol} and a new
##' column in the feature data slot called \code{test} which has part
##' of the labels relabelled as "unknown" class (the number of
##' proteins renamed as "unknown" is according to the parameter size).
##' @seealso \code{\link{sampleMSnSet}} \code{\link{unknownMSnSet}}
##' \code{\link{markerMSnSet}}
##' @author Lisa Breckels
##' @examples
##' library(pRolocdata)
##' data(tan2009r1)
##' sample <- testMSnSet(tan2009r1)
##' getMarkers(sample, "test")
##' all(dim(sample) == dim(markerMSnSet(tan2009r1)))
testMSnSet <- function(object, fcol = "markers",
                       size = .2, seed) {
    if (!missing(seed)) {
        seed <- as.integer(seed)
        set.seed(seed)
    }
    P <- markerMSnSet(object, fcol)
    data <- subsetAsDataFrame(P, fcol, keepColNames = TRUE)
    ## Select validation set
    .size <- ceiling(table(data[ ,fcol]) * size)
    .size <- .size[unique(data[ ,fcol])]
    validation.idxP <- strata(data, fcol, size = .size,
                              method = "srswor")$ID_unit
    validation.names <- rownames(data)[validation.idxP]
    validation.P <- P[validation.names, ]
    train.P <- P[-validation.idxP, ]
    fData(train.P)$test <- as.character(fData(train.P)[, fcol])
    fData(validation.P)$test <- rep("unknown", nrow(validation.P))
    allP <- combine(train.P, validation.P)
    return(allP)
}


##' This function extracts a stratified sample of an \code{MSnSet}.
##' 
##' @title Extract a stratified sample of an \code{MSnSet}
##' @param object An instance of class \code{\linkS4class{MSnSet}}
##' @param fcol The feature meta-data column name containing the
##' marker (vector or matrix) definitions on which the MSnSet will be
##' stratified. Default is \code{markers}.
##' @param size The size of the stratified sample to be
##' extracted. Default is 0.2 (20 percent).
##' @param seed The optional random number generator seed.
##' @return A stratified sample (according to the defined \code{fcol})
##' which is an instance of class \code{"\linkS4class{MSnSet}"}.
##' @seealso \code{\link{testMSnSet}} \code{\link{unknownMSnSet}}
##' \code{\link{markerMSnSet}}. See \code{\link{markers}} for details
##' about markers encoding.
##' @author Lisa Breckels
##' @examples
##' library(pRolocdata)
##' data(tan2009r1)
##' dim(tan2009r1)
##' smp <- sampleMSnSet(tan2009r1, fcol = "markers")
##' dim(smp)
##' getMarkers(tan2009r1)
##' getMarkers(smp)
sampleMSnSet <- function(object, fcol = "markers", size = .2, seed) {
    ## Set seed
    if (!missing(seed)) {
        seed <- as.integer(seed)
        set.seed(seed)
    }
    switch(mrkEncoding(object, fcol),
           vector = vecSampleMSnSet(object, fcol, size),
           matrix = matSampleMSnSet(object, fcol, size))
}

vecSampleMSnSet <- function(object, fcol, size) {
    nms <- sampleNames(object)
    mydata <- data.frame(exprs(object), markers = fData(object)[, fcol])
    colnames(mydata) <- c(nms, fcol)
    subset <- ceiling(table(mydata[, fcol]) * size)
    subset <- subset[unique(mydata[, fcol])]
    idx <- strata(mydata, fcol, size = subset,
                  method = "srswor")$ID_unit
    object <- object[idx,]
    m <- as.character(fData(object)[, fcol])
    tm <- table(m)
    if (any(tm < 6))
        warning("New sample contains classes with < 6 markers",
                call. = FALSE)
    return(object)
}

matSampleMSnSet <- function(object, fcol, size) {
    ## we create a vector of concatenated colnames
    vfcol <- rep("unknown", nrow(object))
    mrk <- fData(object)[, fcol]
    nms <- colnames(mrk)
    for (i in 1:length(vfcol)) {
        k <- nms[mrk[i, ] != 0]
        if (length(k))
            vfcol[i] <- paste(k, collapse = ".")
    }
    ## subsetting
    mydata <- data.frame(exprs(object), markers = vfcol)
    colnames(mydata) <- c(sampleNames(object), fcol)
    subset <- ceiling(table(mydata[, fcol]) * size)
    subset <- subset[unique(mydata[, fcol])]
    idx <- strata(mydata, fcol, size = subset,
                  method = "srswor")$ID_unit
    object <- object[idx, ]
    tm <- table(vfcol)
    if (any(tm < 6))
        warning("New sample contains classes with < 6 markers",
                call. = FALSE)
    return(object)
}

##' Convenience accessor to the organelle classes in an 'MSnSet'.
##' This function returns the organelle classes of an
##' \code{MSnSet} instance. As a side effect, it prints out the classes.
##' 
##' @title Returns the organelle classes in an 'MSnSet'
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The name of the markers column in the \code{featureData}
##' slot. Default is \code{markers}.
##' @param verbose If \code{TRUE}, a character vector of the organelle
##' classes is printed and the classes are returned invisibly. If \code{FALSE}, 
##' the markers are returned.
##' @param ... Additional parameters passed to \code{sort} from the base package.
##' @return A \code{character} vector of the organelle classes in the data.
##' @author Lisa Breckels and Laurent Gatto
##' @seealso \code{\link{getMarkers}} to extract the marker
##' proteins. See \code{\link{markers}} for details about spatial
##' markers storage and encoding.
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' organelles <- getMarkerClasses(dunkley2006)
##' ## same if markers encoded as a matrix
##' dunkley2006 <- mrkVecToMat(dunkley2006, mfcol = "Markers")
##' organelles2 <- getMarkerClasses(dunkley2006, fcol = "Markers")
##' stopifnot(all.equal(organelles, organelles2))
getMarkerClasses <- function(object,
                             fcol = "markers",
                             verbose = TRUE,
                             ...) {
    if (isMrkVec(object, fcol))
        getMarkerVecClasses(object, fcol, verbose, ...)
    else if (isMrkMat(object, fcol))
        getMarkerMatClasses(object, fcol, verbose)
    else
        stop("Your markers are neither vector nor matrix. See ?markers for details.")
}

getMarkerMatClasses <- function(object, fcol, verbose, ...) {
    classes <- colnames(fData(object)[, fcol])
    classes <- sort(classes, ...)
    classes <- classes[which(classes != "unknown")]
    if (verbose) {
        print(classes)
        invisible(classes)
    } else {
        return(classes)
    }
}

getMarkerVecClasses <- function(object, fcol, verbose, ...) {
    organelleMarkers <- getMarkers(object, fcol, verbose = FALSE)
    classes <- unique(organelleMarkers)
    classes <- sort(classes, ...)
    classes <- classes[which(classes != "unknown")]
    if (verbose) {
        print(classes)
        invisible(classes)
    } else {
        return(classes)
    }
}


##' The function assumes that its input is a binary \code{MSnSet} and
##' computes, for each marker class, the number of non-zero expression
##' profiles. The function is meant to be used to produce heatmaps
##' (see the example) and visualise binary (such as GO) \code{MSnSet}
##' objects and assess their utility: all zero features/classes will
##' not be informative at all (and can be filtered out with
##' \code{\link{filterBinMSnSet}}) while features/classes with many
##' annotations (GO terms) are likely not be be informative either. 
##'
##' @title Compute the number of non-zero values in each marker classes
##' @param object An instance of class \code{MSnSet} with binary data.
##' @param fcol A \code{character} defining the feature data variable
##'     to be used as markers. Default is \code{"markers"}.
##' @param as.matrix If \code{TRUE} (default) the data is formatted
##'     and returned as a \code{matrix}. Otherwise, a \code{list} is
##'     returned.
##' @param percent If \code{TRUE}, percentages are
##'     returned. Otherwise, absolute values.
##' @return A \code{matrix} or a \code{list} indicating the number of
##'     non-zero value per marker class.
##' @author Laurent Gatto
##' @seealso \code{\link{filterBinMSnSet}}
##' @examples
##' library(pRolocdata)
##' data(hyperLOPIT2015goCC)
##' zerosInBinMSnSet(hyperLOPIT2015goCC)
##' zerosInBinMSnSet(hyperLOPIT2015goCC, percent = FALSE)
##' pal <- colorRampPalette(c("white", "blue"))
##' library(lattice)
##' levelplot(zerosInBinMSnSet(hyperLOPIT2015goCC),
##'           xlab = "Number of non-0s",
##'           ylab = "Marker class",
##'           col.regions = pal(140))
zerosInBinMSnSet <- function(object, fcol = "markers",
                            as.matrix = TRUE,
                            percent = TRUE) {
    if (!isBinary(object))
        warning("Your assay data is not binary!")
    object <- markerMSnSet(object, fcol = fcol)
    mm <- getMarkerClasses(object, fcol = fcol, verbose = FALSE)
    res <- vector("list", length = length(mm))
    names(res) <- mm
    for (m in mm) {
        mobj <- object[fData(object)[, fcol] == m, ]
        .res <- table(rowSums(exprs(mobj)))
        if (percent) .res <- .res / ncol(mobj)
        res[[m]] <- .res
    }
    if (as.matrix) {
        mres <- matrix(0, ncol = length(res),
                       nrow = max(sapply(res, length)))
        colnames(mres) <- names(res)
        rownames(mres) <- seq_len(nrow(mres)) - 1
        for (m in mm) {
            x <- res[[m]]
            mres[1:length(x), m] <- x
        }
        res <- mres
    }
    res
}

##' Removes columns or rows that have a certain proportion or absolute
##' number of 0 values.
##'
##' @title Filter a binary MSnSet
##' @param object An \code{MSnSet}
##' @param MARGIN 1 or 2. Default is 2.
##' @param t Rows/columns that have \code{t} or less \code{1}s, it
##'     will be filtered out. When \code{t} and \code{q} are missing,
##'     default is to use \code{t = 1}.
##' @param q If a row has a higher quantile than defined by \code{q},
##'     it will be filtered out.
##' @param verbose A \code{logical} defining of a message is to be
##'     printed. Default is \code{TRUE}.
##' @return A filtered \code{MSnSet}.
##' @seealso \code{\link{zerosInBinMSnSet}},
##'     \code{\link{filterZeroCols}}, \code{\link{filterZeroRows}}.
##' @author Laurent Gatto
##' @examples
##' set.seed(1)
##' m <- matrix(sample(0:1, 25, replace=TRUE), 5)
##' m[1, ] <- 0
##' m[, 1] <- 0
##' rownames(m) <- colnames(m) <- letters[1:5]
##' fd <- data.frame(row.names = letters[1:5])
##' x <- MSnSet(exprs = m, fData = fd, pData = fd)
##' exprs(x)
##' ## Remove columns with no 1s
##' exprs(filterBinMSnSet(x, MARGIN = 2, t = 0))
##' ## Remove columns with one 1 or less
##' exprs(filterBinMSnSet(x, MARGIN = 2, t = 1))
##' ## Remove columns with two 1s or less
##' exprs(filterBinMSnSet(x, MARGIN = 2, t = 2))
##' ## Remove columns with three 1s 
##' exprs(filterBinMSnSet(x, MARGIN = 2, t = 3))
##' ## Remove columns that have half or less of 1s
##' exprs(filterBinMSnSet(x, MARGIN = 2, q = 0.5))
filterBinMSnSet <- function(object,
                            MARGIN = 2,
                            t, q,
                            verbose = TRUE) {
    if (!isBinary(object))
        warning("Your assay data is not binary!")
    stopifnot(MARGIN %in% 1:2)
    if (MARGIN == 2)
        K <- colSums(exprs(object))
    else K <- rowSums(exprs(object))
    if (missing(t) & missing(q))
        t <- 1
    if (missing(q)) {
        sel <- K > t
    } else {
        sel <- K > quantile(K, q)
    }
    if (MARGIN == 2) {
        if (verbose) message("Removing ", sum(!sel), " column(s)")
        ans <- object[, sel]
    } else {
        if (verbose) message("Removing ", sum(!sel), " row(s)")
        ans <- object[sel, ]
    }
    if (validObject(ans))
        return(ans)
}

##' Removes all assay data columns/rows that are composed of only 0,
##' i.e. have a \code{colSum}/\code{rowSum} of 0.
##'
##' @title Remove 0 columns/rows
##' @param object A \code{MSnSet} object.
##' @param verbose Print a message with the number of filtered out
##' columns/row (if any).
##' @return An \code{MSnSet}.
##' @author Laurent Gatto

##' @examples
##' library("pRolocdata")
##' data(andy2011goCC)
##' any(colSums(exprs(andy2011goCC)) == 0)
##' exprs(andy2011goCC)[, 1:5] <- 0
##' ncol(andy2011goCC)
##' ncol(filterZeroCols(andy2011goCC))
filterZeroCols <- function(object,
                           verbose = TRUE) {
    cs <- colSums(exprs(object))
    sel <- cs > 0
    if (any(!sel)) {
        if (verbose)
            message("Removing ", sum(!sel), " columns with only 0s.")
        object <- object[, sel]
    }
    if (validObject(object))
        return(object)
}

##' @rdname filterZeroCols
filterZeroRows <- function(object,
                           verbose = TRUE) {
    rs <- rowSums(exprs(object))
    sel <- rs > 0
    if (any(!sel)) {
        if (verbose)
            message("Removing ", sum(!sel), " columns with only 0s.")
        object <- object[sel, ]
    }
    if (validObject(object))
        return(object)
}

##' This function produces organelle-specific quantiles corresponding to 
##' the given classification scores. 
##' 
##' @title Returns organelle-specific quantile scores
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The name of the prediction column in the
##' \code{featureData} slot. 
##' @param scol The name of the prediction score column in the
##' \code{featureData} slot. If missing, created by pasting
##' '.scores' after \code{fcol}. 
##' @param mcol The name of the column containing the training data in the
##' \code{featureData} slot. Default is \code{markers}.
##' @param t The quantile threshold. 
##' @param verbose If \code{TRUE}, the calculated threholds are printed.
##' @return A named \code{vector} of organelle thresholds. 
##' @author Lisa Breckels
##' @seealso \code{\link{getPredictions}} to get organelle predictions based
##' on calculated thresholds.
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' res <- svmClassification(dunkley2006, fcol = "pd.markers",
##'                          sigma = 0.1, cost = 0.5)
##' ## 50% top predictions per class
##' ts <- orgQuants(res, fcol = "svm", t = .5)
##' getPredictions(res, fcol = "svm", t = ts)
orgQuants <- function(object, fcol, scol, 
                      mcol = "markers", 
                      t, verbose = TRUE) {
  stopifnot(!missing(fcol))
  if (missing(scol)) 
    scol <- paste0(fcol, ".scores")
  object <- unknownMSnSet(object, mcol)
  nt <- tapply(fData(object)[, scol], fData(object)[, fcol], quantile, t)
  if (verbose)
    print(nt)
  invisible(nt)
}
