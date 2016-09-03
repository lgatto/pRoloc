setMethod("plot3D", "MSnSet",
          function(object, fcol = "markers", dims = c(1, 2, 3),
                   radius1 = 0.1, radius2 = radius1 * 2, plot = TRUE,
                   ...) {
              x12 <- plot2D(object, dims = dims[1:2], plot = FALSE, ...)
              x13 <- plot2D(object, dims = dims[2:3], plot = FALSE, ...)
              xx <- cbind(x12, x13[, 2, drop = FALSE])
              mm <- factor(getMarkers(object, fcol = fcol, verbose = FALSE))
              ukn <- mm == "unknown"
              col <- getStockcol()
              cls <- col[as.numeric(mm)]
              names(cls) <- mm
              cls[ukn] <- getUnknowncol()
              sel <- !duplicated(names(cls))
              if (plot) {
                  plot3d(xx, col = cls, radius = radius1, type = "s")
                  spheres3d(xx[!ukn, ], radius = radius2, col = cls[!ukn])
                  legend3d("top", legend = names(cls[sel]),
                           col = cls[sel], pch = 19, bty = "n",
                           ncol = 3)
              }
              invisible(xx)
          })

highlightOnPlot3D <- function(object, foi, labels,
                              radius = 0.1 * 3,
                              args = list(), ...) {
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
        text3d(.pca[sel, 1], .pca[sel, 2], .pca[sel, 3], labels[sel], ...)
    } else {
        spheres3d(.pca[sel, 1], .pca[sel, 2], .pca[sel, 3], radius = radius, ...) 
    }
}
