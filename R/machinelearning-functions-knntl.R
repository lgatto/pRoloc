# Sub-functions for `knntlOptimisation` and `knntlClassification`
createPartitions <- function(markers,
                             xval,
                             times,
                             test.size) {
    validation.idx <- caret::createDataPartition(markers, times = times, p = test.size)
    validation <- lapply(validation.idx, function(z) names(markers[z]))
    names(validation) <- paste("validation", 1:times, sep = "")
    train.markers <- lapply(validation.idx, function(z) markers[-z])
    train <- lapply(validation.idx, function(z) names(markers[-z]))
    names(train) <- paste("train", 1:times, sep = "")
    xfolds.idx <- lapply(train.markers, function(z)
        createFolds(z, xval, returnTrain = FALSE))
    train.xfolds <- lapply(1:times,
                           function(x) lapply(xfolds.idx[[x]], function(z) train[[x]][z]))
    names(train.xfolds) <- paste("train", 1:times, sep = "")
    ans <- list(validation = validation,
                train = train,
                train.xfolds = train.xfolds)

    ## Check partitions have been created properly
    .chkPartition <- sapply(1:times, function(z)
        length(unique(union(ans$train[[z]], ans$validation[[z]])))==length(markers))
    if(!all(.chkPartition)) {stop("Validation U training partitions != markers")}
    .chkFolds <- sapply(ans$train.xfolds, function(z) length(unique(unlist(z))))
    .chkFolds <- sapply(1:times, function(z) .chkFolds[z]==length(ans$train[[z]]))
    if(!all(.chkFolds)) {stop("Xfolds are not generated correctly")}

    return(ans)
}


##' The possible weights to be considered is a sequence from 0 (favour
##' auxiliary data) to 1 (favour primary data). Each possible
##' combination of weights for \code{nclass} classes must be
##' tested. The \code{thetas} function produces a weight \code{matrix}
##' for \code{nclass} columns (one for each class) with all possible
##' weight combinations (number of rows).
##'
##' @title Draw matrix of thetas to test
##' @param nclass Number of marker classes
##' @param by The increment of the weights. One of \code{1},
##' \code{0.5}, \code{0.25}, \code{2}, \code{0.1} or \code{0.05}.
##' @param length.out The desired length of the weight sequence.
##' @param verbose A \code{logical} indicating if the weight sequences
##' should be printed out. Default is \code{TRUE}.
##' @return A matrix with all possible theta weight combinations.
##' @author Lisa Breckels
##' @examples
##' dim(thetas(4, by = 0.5))
##' dim(thetas(4, by = 0.2))
##' dim(thetas(5, by = 0.2))
##' dim(thetas(5, length.out = 5))
##' dim(thetas(6, by = 0.2))
thetas <- function(nclass,
                   by = .5,
                   length.out,
                   verbose = TRUE) {
    if (missing(length.out)) {
        bys <- c(1, 0.5, 0.25, 0.2, 0.1, 0.05)
        if (!by %in% bys)
            stop("'by' must be one of ",
                 paste(bys, collapse = ", "))
        t <- seq(0, 1, by)
    } else {
        t <- seq(0, 1, length.out = length.out)
    }
    if (verbose)
        message("Weigths:\n  (", paste(t, collapse = ", "), ")")
    gtools::permutations(length(t), nclass, t, repeats.allowed=TRUE)
}

## getNN function to get va and vp matrices
## include.unknowns = TRUE/FALSE to include unknowns in neighbours
getNN <- function(object, query, labels, k) {
    .classes <- levels(factor(labels))
    if (missing(query)) {
        all <- nrow(object) - 1
        ## Creat matrix of ALL distances
        .nn <- nndist(object, k = all)
    } else {
        all <- nrow(object)
        .nn <- nndist(object, query, k = all)
    }
    .indD <- grep(pattern="dist", x = colnames(.nn), ignore.case = FALSE)
    .indNN <- grep(pattern="index", x = colnames(.nn), ignore.case = FALSE)
    .nnDist <- .nn[, .indD]
    .nnIndex <- .nn[, .indNN]
    nn <- vector("list", nrow(.nnDist))
    for (n in 1:nrow(.nnDist)) {
        if (.nnDist[n, k]==.nnDist[n, k + 1]) {
            tf <- which(.nnDist[n, c((k + 1):all)] != .nnDist[n, k])
            index <- tf[1] + k
            look <- 1:(index - 1)
            nnID <- .nnIndex[n, look]
            nn[[n]] <- labels[nnID]
        } else {
            nnID <- .nnIndex[n, 1:k]
            nn[[n]] <- labels[nnID]
        }
    }
    v <- sapply(.classes, function(x) sapply(nn, function(z) length(which(z==x))))
    rownames(v) <- rownames(.nnDist)
    v <- v/rowSums(v)
    return(v)
}

## Calculates an *individual* vc results given an input known value of theta
vc.res <- function(vp, va, theta) {
    if (length(theta) == ncol(va) & length(theta) == ncol(vp) & ncol(vp) == ncol(va)) {
        n <- 1:ncol(vp)
        .foo <- function(x, y, n) theta[n]*(x) + (1-theta[n])*(y)
        .vc <- sapply(n, function(z) mapply(.foo, x=vp[, z], y=va[, z], n=z))
        colnames(.vc) <- colnames(vp)
        rownames(.vc) <- rownames(vp)
        ## res <- apply(.vc, 1, function(z) names(which.max(z)))
        return(.vc)
    } else {
        stop("Classes differ in input matrices")
    }
}

## Assign class to instance with highest score, if scores are tied pick one
## at random
getPrediction <- function(matrix, seed) {
    if (!missing(seed)) {
        seed <- as.integer(seed)
        set.seed(seed)
    }
    res <- vector("character", length = nrow(matrix))
    for (i in 1:nrow(matrix)) {
        best <- names(which(matrix[i, ] == max(matrix[i, ])))
        if (length(best) > 1) {
            res[i] <- sample(best, size = 1)
        } else {
            res[i] <- best
        }
    }
    return(res)
}

## Function to combine results from splitting theta matrix
combineParams <- function(object) {
    if (!inherits(object, "list")) stop("Object must be a list")

    ## Update algorithm, hyperparams, test.size, design and xval slots
    .algorithm <- object[[1]]@algorithm
    .k <- object[[1]]@hyperparameters$k
    .theta <- lapply(object, function(z) z@hyperparameters$theta)
    .theta <- do.call(rbind, .theta)
    .hyperparams <- list(k = .k,
                         theta = .theta)
    .times <- object[[1]]@design[3] #
    .test.size <- object[[1]]@design[2] #
    .design <- object[[1]]@design
    .xval <- object[[1]]@design[1] #
    .predictions <- .cm <- vector("list", .times)

    ## Combine results and otherWeights over theta splits
    .results <- lapply(object, function(z) z@results)
    .results <- lapply(1:.times,
                       function(x)
                           do.call(rbind, lapply(.results, function(z) z[x, ])))
    .otherWeights <- lapply(object, function(z) z@otherWeights)
    .otherWeights <- lapply(1:.times,
                            function(x)
                                do.call(rbind, lapply(.otherWeights, function(z) z[[x]])))

    ## Get best theta over thetaSplits, if any multiple best theta
    ## across the the thetaSplits rbind these to otherWeights
    for (i in 1:.times) {
        .indMax <- which(.results[[i]][, 1] == max(.results[[i]][, 1]))
        if (length(.indMax) > 1) {
            .choice <- sample(.indMax, 1)
            .otherWeights[[i]] <- rbind(.otherWeights[[i]], .results[[i]][-.choice, -1])
            .indMax <- .choice
        }
        .results[[i]] <- .results[[i]][.indMax, ]
        ## Now get assoicated CM matrix and prediction for best theta over splits
        .predictions[[i]] <- object[[.indMax]]@predictions[[i]]
        .cm[[i]] <- object[[.indMax]]@cmMatrices[[i]]
    }
    .results <- do.call(rbind, .results)

    ## Average f1 scores for each theta
    .f1Matrices <- lapply(1:.times, function(x)
        do.call(c, lapply(object, function(z) z@f1Matrices[[x]])))
    .names <- paste("times", 1:.times, sep="")
    names(.f1Matrices) <- .names
    .thNames <- paste("th", 1:nrow(.theta), sep = "")
    for (i in 1:.times) {names(.f1Matrices[[i]]) <- .thNames}

    ## Save folds and data size
    .folds <- object[[1]]@testPartitions
    .ds <- object[[1]]@datasize

    ## Update log
    if (length(object@log) > 0) {
        .log <- vector("list", 2)
        names(.log) <- c("warnings", "splits")
        .log$warnings <- unlist(lapply(object, function(z) z@log))
        .thList <- sapply(object, function(z) nrow(z@hyperparameters$theta))
        .log$splits <- paste("Theta matrix was split into ", length(object),
                             "matrices for optimal parallelisation with nrow = ",
                             .thList, "for split", 1:length(object))
        names(.log[[1]]) <- NULL
    } else .log <- list()

    ans <- new("ThetaRegRes",
               algorithm = .algorithm,
               hyperparameters = .hyperparams,
               design = .design,
               log = .log,
               results = .results,
               f1Matrices = .f1Matrices,
               cmMatrices = .cm,
               testPartitions = .folds,
               datasize = .ds,
               predictions = .predictions,
               otherWeights = .otherWeights)
    return(ans)
}


## Internal knntlClassification - modified version, with checks removed,
## no generation of MSnSet output
classify <- function(primary,
                     auxiliary,
                     markers,
                     bestTheta,
                     k) {
    ## Get k's
    if (missing(k)) {
        stop("No k passed to function classify")
    } else {
        if(!is.numeric(k)) stop("Input k is not of class 'numeric'")
        if(!length(k)==2) stop("Input k must be of length 2")
    }

    ## Generate nearest neighbours for each protein in primary
    x <- names(markers[which(markers == "unknown")])
    l <- names(markers[which(markers != "unknown")])
    mP <- primary[l, ]
    mA <- auxiliary[l, ]
    uP <- primary[x, ]
    uA <- auxiliary[x, ]
    labels <- markers[l]
    vp <- getNN(mP, uP, labels, k=k[1])
    va <- getNN(mA, uA, labels, k=k[2])

    ## Now select proteins common in both sets
    pn <- rownames(vp)
    an <- rownames(va)
    cmn <- intersect(pn, an)
    pn.idx <- match(cmn, pn)
    an.idx <- match(cmn, an)
    va <- va[an.idx, ]
    vp <- vp[pn.idx, ]
    ## Get vc matrix to vote over
    vcMat <- vc.res(vp, va, bestTheta)
    res <- getPrediction(vcMat)
    names(res) <- rownames(va)
    return(res)
}

## Now split theta matrix
## Generate thetas to use as input
BUG_splitTh <- function(theta, cores) {
    if (cores == 1)
        return(list(theta))
    spl <- round(nrow(theta)/cores)
    t <- vector("list", cores)
    foo <- function(x,y) c(((x-1)*y+1):(x*y))
    for (i in 1:cores) {
        if (i==cores) {
            f <- t[[i-1]][spl]
            f <- f+1
            t[[i]] <- c(f:nrow(theta))
        } else {
            t[[i]] <- foo(i, spl)
        }
    }
    lapply(t, function(z) theta[z, ])
}


splitTh <- function(theta, cores) {
    .checkSplitTh <- function(idxl, n) {
        ## idxl: list with indices
        ## n: nrow(theta)
        tmp <- unlist(idxl)
        all(sort(tmp) == seq_len(n))
    }
    n <- nrow(theta)
    ll <- split(1:n, ceiling(seq_len(n)/(n/cores)))
    stopifnot(.checkSplitTh(ll, n))
    lapply(ll, function(z) theta[z, ])
}


## Core optimisation function for knntlOptimisation.
tlopt <- function(primary,   # matrix
                  auxiliary, # matrix
                  markers,
                  xval,
                  times,
                  k,
                  theta,
                  xfolds) {
    ## set warnings
    .warnings <- NULL
    if (class(theta) == "numeric")
        theta <- t(as.matrix(theta))

    ## get markers, classes and initialise objects
    .f1Matrices <- vector("list", length = times)
    f1Res <- vector("list", xval)

    ## ---------Outer loop for multiple rounds of xval
    for (.times in 1:times) {

        train <- xfolds$train[[.times]]
        validation <- xfolds$validation[[.times]]
        folds <- xfolds$train.xfolds[[.times]]

        ## ----------Inner loop for xval
        for (.xval in 1:xval) {
            .foldNames <- folds[[.xval]]
            .m <- markers[.foldNames]

            trainP <- primary[train, ]
            trainA <- auxiliary[train, ]

            trainMrk <- as.character(markers)
            names(trainMrk) <- names(markers)
            trainMrk <- trainMrk[train]
            trainMrk[.foldNames] <- "unknown"

            llmrk <- levels(as.factor(trainMrk))
            llmrk <- llmrk[which(llmrk != "unknown")]

            ## Classify, get results, calculate confusion matrices and macroF1 scores
            f1Res[[.xval]] <- sapply(1:nrow(theta), function(z) {
                .r <- classify(primary = trainP,
                               auxiliary = trainA,
                               markers = trainMrk,
                               bestTheta = theta[z, ],
                               k = k)
                ## NB: output of classify is only classified unknowns no
                ## markers are included in the output
                .r <- factor(.r, levels = llmrk, ordered = TRUE)
                cm <- table(.r, .m, dnn = c("Prediction", "Reference"))
                f1 <- MLInterfaces:::.macroF1(MLInterfaces:::.precision(cm, naAs0 = TRUE),
                                              MLInterfaces:::.recall(cm, naAs0 = TRUE),
                                              naAs0. = TRUE)})
        }  # -------End inner xval loop
        ## Store f1 means over each round
        .mat <- do.call(rbind, f1Res)
        .f1Matrices[[.times]] <- apply(.mat, 2, function(z) mean(z)) ## Mean of f1s
    }
    return(.f1Matrices)
}


## Function to favour primary data when there is multiple best theta weights
favourPrimary <- function(primary, auxiliary, object,
                          verbose = TRUE) {

    stopifnot(inherits(primary, "MSnSet"))
    stopifnot(inherits(auxiliary, "MSnSet"))
    stopifnot(inherits(object, "ThetaRegRes"))

    .other <- sapply(object@otherWeights, rbind)
    .chk <- sapply(.other, nrow)
    if (all(.chk == 0)) {
        warning("No otherWeights in ThetaRegRes, no possible other best weights")
        return(object)
    }

    fcol <- object@datasize$fcol
    N <- object@design[3]
    k <- object@hyperparameters$k

    ## Initialise new results matrix, otherWeights, confusion matrices
    ncl <- length(getMarkerClasses(primary, fcol))
    results <- matrix(NA, nrow = N, ncol = ncl + 1)
    colnames(results) <- colnames(object@results)

    cmMatrices <- otherWeights <- vector("list", N)

    if (verbose) {
        pb <- txtProgressBar(min = 0,
                             max = N,
                             style = 3)
        ._k <- 0
    }

    for (i in 1:N) {

        if (verbose) {
            setTxtProgressBar(pb, ._k)
            ._k <- ._k + 1
        }

        ow <- object@otherWeights[[i]]
        if (!is.null(ow)) {
            allbest <- rbind(ow, object@results[i, -1])
            rs <- rowSums(allbest)
            indbest <- which(rs == max(rs))

            ## If more than one max(rowSums) sample one
            if (length(indbest) > 1) {
                indbest <- sample(indbest, size = 1)
            }

            th <- allbest[indbest, ]
            otherWeights[[i]] <- allbest[-indbest, ]

            ## Apply best to validation to get macroF1
            val <- object@testPartitions$validation[[i]]
            train <- object@testPartitions$train[[i]]

            primary <- markerMSnSet(primary, fcol)
            auxiliary <- markerMSnSet(auxiliary, fcol)
            auxiliary <- filterZeroCols(auxiliary, verbose = TRUE)

            fData(primary)$xxx <-
                             fData(auxiliary)$xxx <-
                                                as.character(fData(primary)[, fcol])
            fData(primary)[val, "xxx"] <-
                fData(auxiliary)[val, "xxx"] <-
                rep("unknown", length(val))
            res <- knntlClassification(primary, auxiliary, fcol = "xxx",
                                       bestTheta = th, k = k)
            lev <- names(th)

            ## Now calculate macroF1 on validation
            res.x <- unknownMSnSet(res, "xxx")
            r <- factor(getMarkers(res.x, "knntl"),
                        levels = lev)
            m <- factor(getMarkers(res.x, fcol, verbose = FALSE), levels = lev)
            cm <- table(r, m, dnn = c("Prediction", "Reference"))
            f1 <- MLInterfaces:::.macroF1(MLInterfaces:::.precision(cm, naAs0 = TRUE),
                                          MLInterfaces:::.recall(cm, naAs0 = TRUE),
                                          naAs0. = TRUE)

            results[i, ]<- c(f1, th)
            cmMatrices[[i]] <- cm
        } else {
            results[i, ] <- object@results[i, ]
            cmMatrices[[i]] <- object@cmMatrices[[i]]
        }
    }

    if (verbose) {
        setTxtProgressBar(pb, ._k)
        close(pb)
    }

    object@otherWeights <- otherWeights
    object@results <- results
    object@cmMatrices <- cmMatrices

    return(object)
}


##' theta parameter optimisation
##'
##' Classification parameter optimisation for the KNN implementation
##' of Wu and Dietterich's transfer learning schema
##'
##' \code{knntlOptimisation} implements a variation of Wu and
##' Dietterich's transfer learning schema: P. Wu and
##' T. G. Dietterich. Improving SVM accuracy by training on auxiliary
##' data sources. In Proceedings of the Twenty-First International
##' Conference on Machine Learning, pages 871 - 878.  Morgan Kaufmann,
##' 2004. A grid search for the best theta is performed.
##'
##'
##' @param primary An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param auxiliary An instance of class
##'     \code{"\linkS4class{MSnSet}"}.
##' @param fcol The feature meta-data containing marker definitions.
##'     Default is \code{markers}.
##' @param k Numeric vector of length 2, containing the best \code{k}
##'     parameters to use for the primary (\code{k[1]}) and auxiliary
##'     (\code{k[2]}) datasets. See \code{knnOptimisation} for
##'     generating best \code{k}.
##' @param times The number of times cross-validation is
##'     performed. Default is 50.
##' @param test.size The size of test (validation) data. Default is
##'     0.2 (20 percent).
##' @param xval The number of rounds of cross-validation to perform.
##' @param by The increment for theta, must be one of \code{c(1, 0.5,
##'     0.25, 0.2, 0.15, 0.1, 0.05)}
##' @param length.out Alternative to using \code{by}
##'     parameter. Specifies the desired length of the sequence of
##'     theta to test.
##' @param th A matrix of theta values to test for each class as
##'     generated from the function \code{\link{thetas}}, the number
##'     of columns should be equal to the number of classes contained
##'     in \code{fcol}. Note: columns will be ordered according to
##'     \code{getMarkerClasses(primary, fcol)}. This argument is only
##'     valid if the default method 'Breckels' is used.
##' @param xfolds Option to pass specific folds for the cross
##'     validation.
##' @param BPPARAM Required for parallelisation. If not specified
##'     selects a default \code{BiocParallelParam}, from global
##'     options or, if that fails, the most recently registered()
##'     back-end.
##' @param method The k-NN transfer learning method to use. The
##'     default is 'Breckels' as described in the Breckels et al
##'     (2016). If 'Wu' is specificed then the original method
##'     implemented Wu and Dietterich (2004) is implemented.
##' @param log A \code{logical} defining whether logging should be
##'     enabled. Default is \code{FALSE}. Note that logging produes
##'     considerably bigger objects.
##' @param seed The optional random number generator seed.
##' @return A list of containing the theta combinations tested,
##'     associated macro F1 score and accuracy for each combination
##'     over each round (specified by times).
##' @seealso \code{\link{knntlClassification}} and example therein.
##' @aliases knntlOptimisation knntlOptimisation
##' @author Lisa Breckels
##' @references Breckels LM, Holden S, Wonjar D, Mulvey CM,
##'     Christoforou A, Groen AJ, Kohlbacher O, Lilley KS, Gatto L.
##'     Learning from heterogeneous data sources: an application in
##'     spatial proteomics. bioRxiv. doi:
##'     http://dx.doi.org/10.1101/022152
##'
##' Wu P, Dietterich TG. Improving SVM Accuracy by Training on Auxiliary
##' Data Sources. Proceedings of the 21st International Conference on Machine
##' Learning (ICML); 2004.
knntlOptimisation  <- function(primary,
                               auxiliary,
                               fcol = "markers",
                               k,
                               times = 50,
                               test.size = .2,
                               xval = 5,
                               by = .5,
                               length.out,
                               th,
                               xfolds,
                               BPPARAM = BiocParallel::bpparam(),
                               method = "Breckels",
                               log = FALSE,
                               seed) {
    ## Set seed (Originally removed for Darwin HPC, and added
    ## back 22/02/16 to use with SerialParam, unit testing and for
    ## reproducibility)
    if (class(BPPARAM) == "SerialParam") {
        if (missing(seed)) {
            seed <- sample(.Machine$integer.max, 1)
        }
        .seed <- as.integer(seed)
        set.seed(.seed)
    }

    ## Check object validity
    if (!inherits(primary, "MSnSet") | !inherits(auxiliary, "MSnSet"))
        stop("Primary and auxiliary must both be of class 'MSnSet'")

    ## check k specified
    if (missing(k)) {
        stop("No k specified. Generate best k's for primary and auxiliary. See
         ?knnOptimisation")
    } else {
        if (length(k)!=2 | !is.numeric(k)) {
            stop("k must be of class 'numeric' and of length = 2 (one k for each
           data source)")
        }
    }

    ## Check method is valid
    if (!(method == "Breckels" | method == "Wu"))
        stop("method must be one of 'Breckels' or 'Wu'")

    ## We don't care about the unlabelled/unknown instances here, and
    ## want an overlap of marker features. It is not a strict
    ## requirement to have the same markes; different markers will be
    ## used for training, not for validation of testing (to be checked
    ## though).
    primary <- markerMSnSet(primary, fcol)
    auxiliary <- markerMSnSet(auxiliary, fcol)

    ## Filter to remove empty columns - to best tested
    auxiliary <- filterZeroCols(auxiliary, verbose = TRUE)
    classes <- getMarkerClasses(primary, fcol)
    nclass <- length(classes)

    ## From here on don't use MSnSet's, stick to matrices, profiling showed
    ## code was slow using MSnSet's due to class validity checks etc.
    matP <- exprs(primary)
    matA <- exprs(auxiliary)

    mrkP <- as.character(fData(primary)[, fcol])
    mrkA <- as.character(fData(auxiliary)[, fcol])

    ## Check datasets have some common proteins
    if (length(intersect(rownames(matP),
                         rownames(matA))) == 0)
        stop("No common marker proteins in primary and auxilary data")

    if (!identical(sort(unique(mrkP)), sort(unique(mrkA))))
        stop("Different classes in fcol's between data sources")


    ## Generate thetas to test
    if (missing(th)) {
        ## Run 'Wu' method, only data source specific not class-specific
        if (method == "Wu") {
            if (!missing(length.out)) {
                th <- matrix(rep(seq(0, 1, length.out = length.out), nclass),
                             ncol = nclass)
                w1 <- seq(0, 1, length.out)
                w2 <- 1 - w1
            } else {
                th <- matrix(rep(seq(0, 1, by = by), nclass),
                             ncol = nclass)
                w1 <- seq(0, 1, length.out)
                w2 <- 1 - w1
            }
            colnames(th) <- classes
            .numTh <- nrow(th)
        } else {
            ## Run 'Breckels' method, generate all possible theta weights
            if (!missing(length.out)) {
                th <- thetas(nclass, length.out = length.out)
                w1 <- seq(0, 1, length.out = length.out)
                w2 <- 1 - w1
            } else {
                th <- thetas(nclass, by)
                w1 <- seq(0, 1, by)
                w2 <- 1 - w1
            }
            colnames(th) <- classes
            .numTh <- nrow(th)
        }
    } else {
        if (method != "Breckels")
            stop("Wu's method can not be used with a specific theta matrix,
           please run with method = 'Breckels'")
        ## Check the input matrix
        nP <- names(table(fData(primary)[,fcol]))
        nA <- names(table(fData(auxiliary)[,fcol]))
        if (any(nP == "unknown")) {
            nP <- length(nP)-1
            nA <- length(nA)-1
        } else {
            nP <- length(nP)
            nA <- length(nA)
        }
        if (!is.matrix(th)) stop("thetas is not a matrix")

        if (!all(nP == ncol(th)))
            stop("thetas has a different number of classes to primary data")
        if (!all(nA == ncol(th)))
            stop("thetas has a different number of classes to primary data")
        if (is.null(colnames(th))) {
            message("Note: vector will be ordered according to classes: ",
                    paste(classes, collpase = ""),
                    "(as names are not explicitly defined)")
            colnames(th) <- classes
        } else {
            th <- th[, c(match(classes, colnames(th))), drop = FALSE]
        }
        .numTh <- nrow(th)
        w1 <- unique(as.vector(th))
        w2 <- 1-w1
    }

    ## Now select proteins common in both sets
    pn <- rownames(matP)
    an <- rownames(matA)
    cmn <- intersect(pn, an)
    pn.idx <- match(cmn, pn)
    an.idx <- match(cmn, an)
    matP <- matP[pn.idx, ]
    matA <- matA[an.idx, ]
    mrkP <- mrkP[pn.idx]
    mrkA <- mrkA[an.idx]
    names(mrkP) <- rownames(matP)
    names(mrkA) <- rownames(mrkA)

    markers <- mrkP

    ## Do subsetting for train/validation partition if not specified in xfolds
    if (missing(xfolds)) {
        xfolds <- createPartitions(markers, xval, times, test.size)
    } else {
        if(class(xfolds) != "list" | length(xfolds) != 3 |
           length(xfolds[[1]]) != times | class(xfolds[[1]][[1]])!="character") {
            stop("If 'xfolds' is specified it must be generated using ",
                 "the function pRoloc::createPartitons")
        }
    }

    .workers <- as.numeric(BPPARAM$workers)
    if (.numTh < .workers) {
        .workers <- .numTh # Is there enough rows in the matrix to split amongst cores
    }
    .thetaSubsets <- splitTh(theta = th, cores = .workers)
    .res <- bplapply(.thetaSubsets,
                     function(z) {
                         tlopt(primary = matP,
                               auxiliary = matA,
                               markers = markers,
                               ## fcol = fcol,
                               xval = xval,
                               times = times,
                               k = k,
                               theta = z,
                               xfolds = xfolds
                               ## test.size = test.size
                               )},
                     BPPARAM = BPPARAM)
    .f1Matrices <- sapply(1:times, function(x)
        unlist(lapply(.res, function(z) z[[x]])), simplify = FALSE)

    ## Initialise objects for ThetaRegRes
    .thNames <- paste("th", 1:nrow(th), sep = "")

    .otherWeights <- .cmMatrices <- .pred <- vector("list", times)
    results <- matrix(NA, nrow = times, ncol = nclass + 1)
    .warnings <- NULL

    for (.times in 1:times) {
        names(.f1Matrices[[.times]]) <- .thNames
        ## Getting the best theta
        if (.numTh != 1) {
            ## Find thetas with highest macroF1
            allbest <- which(.f1Matrices[[.times]] == max(.f1Matrices[[.times]]))

            ## Are there many best thetas?
            if (length(allbest) > 0) {
                if (log)
                    .warnings <- c(.warnings,
                                   paste0("For run number ", .times,
                                          " of times there are multiple best thetas, ",
                                          "picking one at random"))
                indbest <- sample(1:length(allbest), 1)
                indTh <- allbest[indbest]
                indOthers <- allbest[-indbest]
                .otherWeights[[.times]] <- th[indOthers, ]
            } else {
                indTh <- allbest
                ## .otherWeights[[.times]] <- NULL
            }
            .best <- th[indTh, ]
        } else {
            ## a matrix with 1 row
            .best <- as.numeric(th)
            names(.best) <- colnames(th)
        }

        ## Apply best to validation to get macroF1
        val <- xfolds$validation[[.times]]
        train <- xfolds$train[[.times]]

        P <- primary
        A <- auxiliary

        fData(P)$xxx <- as.character(fData(P)[, fcol])
        fData(A)$xxx <- as.character(fData(A)[, fcol])
        fData(P)[val, "xxx"] <- fData(A)[val, "xxx"] <- rep("unknown", length(val))
        res <- knntlClassification(P, A, fcol = "xxx",
                                   bestTheta = .best, k = k)
        lev <- names(.best)

        ## Now calculate macroF1 on validation
        res.x <- unknownMSnSet(res, "xxx")
        .pred[[.times]] <- factor(getMarkers(res.x, "knntl", verbose = FALSE),
                                  levels = lev)
        .mark <- factor(getMarkers(res.x, fcol, verbose = FALSE), levels = lev)
        cm <- table(.pred[[.times]], .mark, dnn = c("Prediction", "Reference"))
        f1 <- MLInterfaces:::.macroF1(MLInterfaces:::.precision(cm, naAs0 = TRUE),
                                      MLInterfaces:::.recall(cm, naAs0 = TRUE),
                                      naAs0. = TRUE)

        results[.times, ]<- c(f1, .best)
        .cmMatrices[[.times]] <- cm
    }
    colnames(results) <- c("F1", lev)

    ## Get datasizes information to pass to ThetaRegRes slot
    trainP1 <- markerMSnSet(P, "xxx")
    .idxNames <- xfolds$train.xfolds[[times]][[1]]
    .idx <- match(.idxNames, featureNames(trainP1))
    trainP2 <- trainP1[-.idx, ]
    trainA1 <- markerMSnSet(A, "xxx")
    .idxA <- match(.idxNames, featureNames(trainA1))
    trainA2 <- trainP1[-.idxA, ]

    ## Store everything else for 'thetaRegRes' instance
    .hyperparams <- list(k = k, theta = th)
    .design <- c(xval = xval, test.size = test.size, times = times)
    .ds <- list(
        "primary" = dim(primary),
        "train.primary1" = dim(trainP1),
        "validation.primary1" = c(length(val), ncol(P)),
        "train.primary2" = dim(trainP2),
        "validation.primary2" = c(length(.idx), ncol(P)),
        "primary.markers" = table(fData(primary)[, fcol]),
        "train.primary.markers1" = table(fData(trainP1)[, fcol]),
        "train.primary.markers2" = table(fData(trainP2)[, fcol]),
        "auxiliary" = dim(auxiliary),
        "train.auxiliary1" = dim(trainA1),
        "validation.auxiliary1" = c(length(val), ncol(A)),
        "train.auxiliary2" = dim(trainA2),
        "validation.auxiliary2" = c(length(.idxA), ncol(A)),
        "auxiliary.markers" = table(fData(auxiliary)[, fcol]),
        "train.auxiliary.markers1" = table(fData(trainA1)[, fcol]),
        "train.auxiliary.markers2" = table(fData(trainA2)[, fcol]),
        "fcol" = fcol
    )
    .xfolds <- xfolds[-3]

    ## Update log
    .log <- list()
    if (log) {
        .log <- vector("list", 2)
        names(.log) <- c("warnings", "splits")
        if (is.null(.warnings)) {
            .log$warnings <- NA
        } else {
            .log$warnings <- .warnings
        }
        .thList <- sapply(.thetaSubsets, nrow)
        .log$splits <- sapply(.thList, function(z)
            paste("Theta matrix was split into ", .workers,
                  "matrices for optimal parallelisation with nrow = ", z))
    }
    if (!log)
        .f1Matrices <- list()
    ans <- new("ThetaRegRes",
               algorithm = "theta",
               hyperparameters = .hyperparams,
               design = .design,
               results = results,
               f1Matrices = .f1Matrices,
               cmMatrices = .cmMatrices,
               testPartitions = .xfolds,
               datasize = .ds,
               otherWeights = .otherWeights,
               predictions = .pred,
               log = .log)
    return(ans)
}

##' Classification using a variation of the KNN implementation
##' of Wu and Dietterich's transfer learning schema
##'
##' @title knn transfer learning classification
##' @param primary An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param auxiliary An instance of class
##'     \code{"\linkS4class{MSnSet}"}.
##' @param fcol The feature meta-data containing marker definitions.
##'     Default is \code{markers}.
##' @param bestTheta Best theta vector as output from
##'     \code{knntlOptimisation}, see \code{knntlOptimisation} for
##'     details
##' @param k Numeric vector of length 2, containing the best \code{k}
##'     parameters to use for the primary and auxiliary datasets. If k
##'     \code{k} is not specified it will be calculated internally.
##' @param scores One of \code{"prediction"}, \code{"all"} or
##'     \code{"none"} to report the score for the predicted class
##'     only, for all classes or none.
##' @param seed The optional random number generator seed.
##' @return A character vector of the classifications for the unknowns
##' @seealso \code{\link{knntlOptimisation}}
##' @author Lisa Breckels
##' @examples
##' \donttest{
##' library(pRolocdata)
##' data(andy2011)
##' data(andy2011goCC)
##' ## reducing calculation time of k by pre-running knnOptimisation
##' x <- c(andy2011, andy2011goCC)
##' k <- lapply(x, function(z)
##'             knnOptimisation(z, times=5,
##'                             fcol = "markers.orig",
##'                             verbose = FALSE))
##' k <- sapply(k, function(z) getParams(z))
##' k
##' ## reducing parameter search with theta = 1,
##' ## weights of only 1 or 0 will be considered
##' opt <- knntlOptimisation(andy2011, andy2011goCC,
##'                          fcol = "markers.orig",
##'                          times = 2,
##'                          by = 1, k = k)
##' opt
##' th <- getParams(opt)
##' plot(opt)
##' res <- knntlClassification(andy2011, andy2011goCC,
##'                            fcol = "markers.orig", th, k)
##' res
##' }
knntlClassification <- function(primary,
                                auxiliary,
                                fcol = "markers",
                                bestTheta,
                                k,
                                scores = c("prediction", "all", "none"),
                                seed) {

    ## Set seed
    if (missing(seed)) {
        seed <- sample(.Machine$integer.max, 1)
    }
    .seed <- as.integer(seed)
    set.seed(.seed)

    scores <- match.arg(scores)
    if (!inherits(primary, "MSnSet") | !inherits(auxiliary, "MSnSet"))
        stop("Primary and auxiliary must both be of class 'MSnSet'")

    markers <- getMarkers(primary, fcol, verbose = FALSE)
    classes <- getMarkerClasses(primary, fcol)
    if (!any(markers == "unknown"))
        stop("No unknown proteins to classify")

    ## In knntlClassification, we want to have the same unknowns [*],
    ## not necessarily in the same order. Different markers are not a
    ## problem (although this should not happen, as not allowed in
    ## knntlOptimisation, where some overlap in needed). We look at
    ## the unknowns' nearest neighbours independently.  [*] it would
    ## be possible to have different ones, but we don't bother.
    if (!checkSortedFeatureNames(unknownMSnSet(primary, fcol),
                                 unknownMSnSet(auxiliary, fcol)))
        stop("Feature names of unknown features don't match exactly.")

    if (inherits(bestTheta, "ThetaRegRes")) {
        k <- bestTheta@hyperparameters$k
        bestTheta <- getParams(bestTheta)
    }

    if (length(bestTheta) != length(classes))
        stop("Classes in best theta and classes in data do not match")

    if (!identical(names(bestTheta), classes))
        stop("Classes in data and bestTheta do not match")

    ## Get k's
    if (missing(k)) {
        stop("Use the same k as for knntlOptimisation.")
    } else {
        if (!is.numeric(k)) stop("Input k is not of class 'numeric'")
        if (!length(k) == 2) stop("Input k must be of length 2")
    }


    ## Generate nearest neighbours for each protein in primary
    matP <- exprs(primary)
    matA <- exprs(auxiliary)
    x <- names(markers[which(markers == "unknown")])
    l <- names(markers[which(markers != "unknown")])
    mP <- matP[l, ]
    mA <- matA[l, ]
    uP <- matP[x, ]
    uA <- matA[x, ]
    labels <- markers[l]
    vp <- getNN(mP, uP, labels, k=k[1])
    va <- getNN(mA, uA, labels, k=k[2])

    ## Now select proteins common in both sets
    pn <- rownames(vp)
    an <- rownames(va)
    cmn <- intersect(pn, an)
    pn.idx <- match(cmn, pn)
    an.idx <- match(cmn, an)
    va <- va[an.idx, ]
    vp <- vp[pn.idx, ]

    ## Get vc matrix to vote over
    vcMat <- vc.res(vp, va, bestTheta)

    ## Match labelled and unlabelled with original MSnSet indices
    L <- match(rownames(mP), featureNames(primary))
    X <- match(rownames(uP), featureNames(primary))

    ## Added as check for LMB (will remove later)
    if (all(colnames(vcMat) != classes))
        stop("Column names in vote matrix not equal to classes")
    if (scores == "all") {
        .scoreMat <- matrix(data = NA, nrow = nrow(primary),
                            ncol = length(classes))
        colnames(.scoreMat) <- paste0(colnames(vcMat), ".knntl.scores")
        for (i in 1:length(classes)) {
            .ind <- which(fData(primary)[, fcol] == classes[i])
            .mark <- rep(1, length(.ind))
            .un <- rep(0, nrow(primary))
            .un[.ind] <- .mark
            .scoreMat[, i] <- .un
        }
        .scoreMat[X, ] <- vcMat
        fData(primary)$knntl.all.scores <- .scoreMat
    }
    else if (scores == "prediction") {
        scores <- apply(vcMat, 1, function(z) max(z))
        knntl.scores <- vector("numeric", nrow(primary))
        knntl.scores[L] <- rep(1, length(L))
        knntl.scores[X] <- scores
        fData(primary)$knntl.scores <- knntl.scores
    }

    ## Get final classification
    res <- getPrediction(vcMat)
    y <- rep("unknown", nrow(primary))
    y[L] <- as.character(fData(primary)[L, fcol])
    y[X] <- res
    fData(primary)$knntl <- as.factor(y)

    return(primary)
}


## Function to empty slots that are not populated during
## knntlOptimisation(..., log = FALSE).
nologgin <- function(x) {
    stopifnot(inherits(x, "ThetaRegRes"))
    x@log <- x@f1Matrices <- list()
    if (validObject(x)) return(x)
}
