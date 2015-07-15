.makeCol3 <- function(x, y, fcol) {
    col0 <- getStockcol()
    setStockcol(NULL)
    col <- paste0(getStockcol(), 90)
    .fcol1 <- factor(fData(x)[, fcol])
    .fcol2 <- factor(fData(y)[, fcol])
    col1 <- col[as.numeric(.fcol1)]
    col2 <- col[as.numeric(.fcol2)]
    col1[.fcol1 == "unknown"] <-
        paste0(getUnknowncol(), 30)
    col2[.fcol2 == "unknown"] <-
        paste0(getUnknowncol(), 30)
    return(list(col0, col1, col2))
}


##' Takes 2 \code{linkS4class{MSnSet}} instances as input to plot the
##' two data sets on the same PCA plot. The second data points are
##' projected on the PC1 and PC2 dimensions calculated for the first
##' data set.
##'
##' @title Draw 2 data sets on one PCA plot
##' @param object An \code{linkS4class{MSnSet}} or a
##' \code{MSnSetList}. In the latter case, only the two first elements
##' of the list will be used for plotting and the others will be
##' silently ignored.
##' @param pcol If \code{object} is an \code{MSnSet}, a \code{factor}
##' or the name of a phenotype variable (\code{phenoData} slot)
##' defining how to split the single \code{MSnSet} into two or more
##' data sets.  Ignored if \code{object} is a
##' \code{\link{MSnSetList}}.
##' @param fcol Feature meta-data label (fData column name) defining
##' the groups to be differentiated using different colours. Default
##' is \code{markers}. Use \code{NULL} to suppress any colouring.
##' @param cex.x Character expansion for the first data set. Default
##' is 1.
##' @param cex.y Character expansion for the second data set. Default
##' is 1.
##' @param pch.x Plotting character for the first data set. Default is
##' 21.
##' @param pch.y Plotting character for the second data set. Default
##' is 23.
##' @param col A vector of colours to highlight the different classes
##' defined by \code{fcol}. If missing (default), default colours are
##' used (see \code{\link{getStockcol}}).
##' @param mirrorX A \code{logical} indicating whether the x axis
##' should be mirrored?
##' @param mirrorY A \code{logical} indicating whether the y axis
##' should be mirrored?
##' @param plot If \code{TRUE} (default), a plot is produced.
##' @param ... Additinal parameters passed to \code{plot} and
##' \code{points}.
##' @return Used for its side effects of producing a plot. Invisibly
##' returns an object of class \code{plot2Ds}, which is a list with
##' the PCA analsys results (see \code{\link{prcomp}}) of the first
##' data set and the new coordinates of the second data sets, as used
##' to produce the plot and the respective point colours. Each of
##' these elements can be accessed with \code{data1}, \code{data2},
##' \code{col1} and \code{code2} respectively.
##' @author Laurent Gatto
##' @seealso See \code{\link{plot2D}} to plot a single data set and
##' \code{\link{move2Ds}} for a animation.
##' @aliases data1 data2 col1 col2
##' @examples
##' library("pRolocdata")
##' data(tan2009r1)
##' data(tan2009r2)
##' msnl <- MSnSetList(list(tan2009r1, tan2009r2))
##' plot2Ds(msnl)
##' ## tweaking the parameters
##' plot2Ds(list(tan2009r1, tan2009r2),
##'         fcol = NULL, cex.x = 1.5)
##' ## input is 1 MSnSet containing 2 data sets
##' data(dunkley2006)
##' plot2Ds(dunkley2006, pcol = "replicate")
##' ## no plot, just the data
##' res <- plot2Ds(dunkley2006, pcol = "replicate",
##'                plot = FALSE)
##' res
##' head(data1(res))
##' head(col1(res))
plot2Ds <- function(object, pcol,
                    fcol = "markers",
                    cex.x = 1, cex.y = 1,
                    pch.x = 21, pch.y = 23,
                    col,
                    mirrorX = FALSE,
                    mirrorY = FALSE,
                    plot = TRUE,
                    ...) {
    if (inherits(object, "MSnSet")) {
        if (missing(pcol)) stop("specify an pcol.")
        object <- split(object, f = pcol)
    } ## a list or an MSnSetList
    x <- object[[1]]
    y <- object[[2]]
    ## prepare the data
    pca1 <- prcomp(exprs(x), scale = TRUE, center = TRUE)
    newdata <- exprs(y)
    colnames(newdata) <- sampleNames(x)
    pred2 <- predict(pca1, newdata = newdata)
    res <- list(pca = pca1, pred = pred2)
    if (missing(col) & !is.null(fcol)) {
        lcols <- .makeCol3(x, y, fcol)
        col0 <- lcols[[1]]
        res$pca.col <- col1 <- lcols[[2]]
        res$pred.col <- col2 <- lcols[[3]]
    }
    class(res) <- c("plot2Ds", "list")
    if (plot) { ## plotting
        if (mirrorX) {
            pca1$x[, 1] <- -pca1$x[, 1]
            pred2[, 1] <- -pred2[, 1]
        }
        if (mirrorY) {
            pca1$x[, 2] <- -pca1$x[, 2]
            pred2[, 2] <- -pred2[, 2]
        }
        xlim <- range(c(pca1$x[, 1], pred2[, 1]))
        ylim <- range(c(pca1$x[, 2], pred2[, 2]))
        ## First plot unknown
        if (is.null(fcol)) {
            unk1 <- unk2 <- TRUE
            col1 <- col2 <- paste0(getUnknowncol(), 30)
        } else {
            unk1 <- fData(x)[, fcol] == "unknown"
            unk2 <- fData(y)[, fcol] == "unknown"
        }
        plot(pca1$x[unk1, 1:2], bg = col1[unk1],
             pch = pch.x,
             xlim = xlim, ylim = ylim,
             cex = cex.x, ...)
        points(pred2[unk2, 1], pred2[unk2, 2],
               pch = pch.y, bg = col2[unk2],
               cex = cex.y)
        grid()
        ## segments
        segments(pca1$x[, 1],
                 pca1$x[, 2],
                 pred2[, 1],
                 pred2[, 2],
                 lty = "dotted",
                 col = "#000000AA")
        ## Add other points
        points(pca1$x[!unk1, 1:2], bg = col1[!unk1],
               pch = pch.x,
               xlim = xlim, ylim = ylim,
               cex = cex.x, ...)
        points(pred2[!unk2, 1], pred2[!unk2, 2],
               pch = pch.y, bg = col2[!unk2],
               cex = cex.y)
        if (!is.null(fcol)) setStockcol(col0)
    }
    invisible(res)
}

print.plot2Ds <- function(x, ...) {
    cat("A 'plot2Ds' result\n")
    cat(" dim(data1): [", nrow(x$pca$x), ", ", ncol(x$pca$x), "]\n", sep = "")
    cat(" dim(data2): [", nrow(x$pred), ", ", ncol(x$pred), "]\n", sep = "")
}

data1 <- function(x) x$pca$x[, 1:2]

data2 <- function(x) x$pred[, 1:2]

col1 <- function(x) x$pca.col

col2 <- function(x) x$pred.col


##' Given two \code{MSnSet} instances of one \code{MSnSetList} with at
##' least two items, this function produces an animation that shows
##' the transition from the first data to the second.
##'
##' @title Displays a spatial proteomics animation
##' @param object An \code{linkS4class{MSnSet}} or a
##' \code{MSnSetList}. In the latter case, only the two first elements
##' of the list will be used for plotting and the others will be
##' silently ignored.
##' @param pcol If \code{object} is an \code{MSnSet}, a \code{factor}
##' or the name of a phenotype variable (\code{phenoData} slot)
##' defining how to split the single \code{MSnSet} into two or more
##' data sets.  Ignored if \code{object} is a
##' \code{\link{MSnSetList}}.
##' @param fcol Feature meta-data label (fData column name) defining
##' the groups to be differentiated using different colours. Default
##' is \code{markers}. Use \code{NULL} to suppress any colouring.
##' @param n Number of frames, Default is 25.
##' @param hl An optional instance of class
##' \code{linkS4class{FeaturesOfInterest}} to track features of
##' interest.
##' @return Used for its side effect of producing a short animation.
##' @seealso \code{\link{plot2Ds}} to a single figure with the two
##' datasets.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' 
##' ## Create a relevant MSnSetList using the dunkley2006 data
##' xx <- split(dunkley2006, "replicate")
##' xx1 <- xx[[1]]
##' xx2 <- xx[[2]]
##' fData(xx1)$markers[374] <- "Golgi"
##' fData(xx2)$markers[412] <- "unknown"
##' xx@@x[[1]] <- xx1
##' xx@@x[[2]] <- xx2
##'
##' ## The features we want to track
##' foi <- FeaturesOfInterest(description = "test",
##'                           fnames = featureNames(xx[[1]])[c(374, 412)])
##'
##' ## (1) visualise each experiment separately
##' par(mfrow = c(2, 1))
##' plot2D(xx[[1]], main = "condition A")
##' highlightOnPlot(xx[[1]], foi)
##' plot2D(xx[[2]], mirrorY = TRUE, main = "condition B")
##' highlightOnPlot(xx[[2]], foi, args = list(mirrorY = TRUE))
##'
##' ## (2) plot both data on the same plot
##' par(mfrow = c(1, 1))
##' tmp <- plot2Ds(xx) 
##' highlightOnPlot(data1(tmp), foi, lwd = 2)
##' highlightOnPlot(data2(tmp), foi, pch = 5, lwd = 2)
##'
##' ## (3) create an animation
##' move2Ds(xx, pcol = "replicate")
##' move2Ds(xx, pcol = "replicate", hl = foi)
move2Ds <- function(object, pcol,
                    fcol = "markers",
                    n = 25,
                    hl) {
    x <- pRoloc::plot2Ds(object, pcol, fcol, plot = FALSE)
    start <- x$pca$x[, 1:2]
    end <- x$pred[, 1:2]
    dx <- end[, 1] - start[, 1]
    dy <- end[, 2] - start[, 2]
    xlim <- range(c(start[, 1], end[, 1]))
    ylim <- range(c(start[, 2], end[, 2]))

    hlidx <- FALSE
    if (!missing(hl))
        hlidx <- rownames(start) %in% foi(hl)

    if (is.null(fcol)) {
        colpar <- matrix(getUnknowncol(), ncol = n, nrow = nrow(start))
    } else {
        colpar <- matrix(NA_character_, ncol = n, nrow = nrow(start))
        for (i in 1:nrow(start))
            colpar[i, ] <- colorRampPalette(c(x$pca.col[i], x$pred.col[i]))(n)
    }
    kk <- 1
    for (k in seq(0, 1, length = n)) {
        X <- Y <- numeric(nrow(start))
        for (i in 1:nrow(start)) {
            X[i] <- start[i, 1] + (k * dx[i])
            Y[i] <- start[i, 2] + (k * dy[i])
        }
        plot(X, Y, xlim = xlim, ylim = ylim,
             xlab = "PC1", ylab = "PC2",
             bg = colpar[, kk], pch = 21)
        segments(start[, 1], start[, 2],
                 end[, 1], end[, 2],
                 col = "#00000020")
        points(X[hlidx], Y[hlidx], cex = 3, lwd = 2)
        kk <- kk + 1
    }
}
