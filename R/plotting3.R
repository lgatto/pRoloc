##' @param radius1 A \code{numeric} specifying the radius of feature
##'     of unknown localisation. Default is 0.1, which is specidied on
##'     the data scale. See \code{\link[rgl]{plot3d}} for details.
##' @param radius2 A \code{numeric} specifying the radius of marker
##'     feature. Default is \code{radius} * 2.
##' @rdname plot2D
##' @examples
##'
##' ## plotting in 3 dimenstions
##' plot3D(dunkley2006)
##' plot3D(dunkley2006, radius2 = 0.3)
##' plot3D(dunkley2006, dims = c(2, 4, 6))
setMethod("plot3D", "MSnSet",
          function(object, fcol = "markers", dims = c(1, 2, 3),
                   radius1 = 0.1, radius2 = radius1 * 2, plot = TRUE,
                   ...) {
              if (!require("rgl"))
                  stop(paste0("Plotting in 3D depends on the rgl package.\n",
                              "Install it with 'install.packages('rgl')'."))
              x12 <- plot2D(object, dims = dims[1:2], plot = FALSE, fcol = fcol, ...)
              x13 <- plot2D(object, dims = dims[2:3], plot = FALSE, fcol = fcol, ...)
              xx <- cbind(x12, x13[, 2, drop = FALSE])
              mm <- factor(getMarkers(object, fcol = fcol, verbose = FALSE))
              ukn <- mm == "unknown"
              col <- getStockcol()
              cls <- col[as.numeric(mm)]
              names(cls) <- mm
              cls[ukn] <- getUnknowncol()
              sel <- !duplicated(names(cls))
              if (plot) {
                  rgl::plot3d(xx, col = cls, radius = radius1, type = "s")
                  rgl::spheres3d(xx[!ukn, ], radius = radius2, col = cls[!ukn])
                  rgl::legend3d("top", legend = names(cls[sel]),
                                col = cls[sel], pch = 19, bty = "n",
                                ncol = 3)
              }
              invisible(xx)
          })
##' @param radius Radius of the spheres to be added to the
##'     visualisation produced by \code{plot3D}. Default is 0.3 (i.e
##'     \code{plot3D}'s \code{radius1} * 3), to emphasise the features
##'     with regard to uknown (\code{radius1 = 0.1}) and marker
##'     (\code{radius1} * 2) features.
##' @rdname highlightOnPlot
##' @examples
##'
##' ## in 3 dimensions
##' if (interactive()) {
##'   plot3D(tan2009r1, radius1 = 0.05)
##'   highlightOnPlot3D(tan2009r1, x, labels = TRUE)
##'   highlightOnPlot3D(tan2009r1, x)
##' }
highlightOnPlot3D <- function(object, foi, labels,
                              args = list(),
                              radius = 0.1 * 3,
                              ...) {
    requireNamespace("rgl")
    if (is.character(foi))
        foi <- FeaturesOfInterest(description = "internally created",
                                  fnames = foi)
    stopifnot(inherits(foi, "FeaturesOfInterest"))
    if (!fnamesIn(foi, object)) {
        warning("None of the features of interest are present in the data.")
        return(invisible(NULL))
    }
    if (inherits(object, "MSnSet")) {
        .args <- list(object = object, plot = FALSE)
        args <- c(args, .args)
        .pca <- do.call(plot3D, args = args)
        sel <- featureNames(object) %in% foi(foi)
    } else if (is.matrix(object)) {
        .pca <- object
        sel <- rownames(object) %in% foi(foi)
    } else {
        stop("'object' must be a matrix (as returned by plot2D) or an MSnSet.")
    }
    if (!missing(labels)) {
        if (is.character(labels)) {
            stopifnot(inherits(object, "MSnSet"))
            labels <- labels[1]
            stopifnot(labels %in% fvarLabels(object))
            labels <- fData(object)[, labels]
        } else if (isTRUE(labels)) {
            if (inherits(object, "MSnSet"))
                labels <- featureNames(object)
            else labels <- rownames(object) ## a matrix
        } else {
            stop("'labels' must be a character or logical of length 1.")
        }
        rgl::text3d(.pca[sel, 1], .pca[sel, 2], .pca[sel, 3], labels[sel], ...)
    } else {
        rgl::spheres3d(.pca[sel, 1], .pca[sel, 2], .pca[sel, 3],
                       radius = radius, ...)
    }
}
