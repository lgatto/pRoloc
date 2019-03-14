## metropolisHastings and plotMetropolisHastings are private functions used to
## generate illustrative MCMC cartoons (see the TAGM workflow as an example).
##
## Usage:
##
## set.seed(2)
## met <- metropolisHastings(10000, 0.9)
##
## par(mfrow = c(2, 2),
##     mar = c(4, 4, 2, 1))
##
## plotMetropolisHastings(2, met)
## plotMetropolisHastings(3, met)
## plotMetropolisHastings(5, met)
## plotMetropolisHastings(30, met)

metropolisHastings <- function (n, rho) {
    mat <- matrix(ncol = 2, nrow = n)   ## matrix for storing the random samples
    x <- y <- 4   ## initial values for all parameters
    cov <- matrix(c(1,sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2)
    prev <- dmvnorm(c(x, y), mean = c(0, 0), sigma = cov)
    mat[1, ] <- c(x, y)  ## initialize the markov chain
    counter <- 1
    while (counter <= n) {
        propx <- x + rnorm(1, 0, 0.1)
        propy <- y + rnorm(1, 0, 0.1)

        newprob <- mixtools::dmvnorm(c(propx, propy), sigma = cov)
        ratio <- newprob/prev

        prob.accept <- min(1, ratio) ## ap
        rand <- runif(1)
        if (rand <= prob.accept) {
            x <- propx
            y <- propy
            mat[counter, ] <- c(x, y) ## store this in the storage array
            counter <- counter + 1
            prev <- newprob ## get ready for the next iteration
        }

    }
    return(mat)
}


plotMetropolisHastings <- function(r, met, rho) {
    mycolb <- rgb(0, 0, 255, maxColorValue = 255, alpha = 125, names = "blue50")
    mycolr <- rgb(255, 0, 0, maxColorValue = 255, alpha = 175, names = "red50")
    x <- mixtools::rmvnorm(10000,
                           sigma = matrix(c(1, sqrt(1 - rho^2), sqrt(1 - rho^2), 1),
                                          ncol = 2))
    a <- x[1:35, 1:2]
    plot(a, ylim = c(-5,5), xlim = c(-5,5),
         xlab = "Channel 1", ylab = "Channel 2",
         col = mycolb, cex = 2, pch = 19)
    legend("topleft", legend = paste0("Iteration ", r), bty = "n")
    mixtools::ellipse(mu = c(0, 0),
                      sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2),
                      alpha = .05, npoints = 1000, newplot = FALSE, draw = TRUE)
    mixtools::ellipse(mu = c(0, 0),
                      sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2),
                      alpha = .01, npoints = 1000, newplot = FALSE, draw = TRUE)
    mixtools::ellipse(mu = c(0, 0),
                      sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2),
                      alpha = .10, npoints = 1000, newplot = FALSE, draw = TRUE)

    points(met[100 * c(1:r) - 99,],
           type = "b", col = mycolr,
           cex = 2, lwd = 3, pch = 19)
}


plotDist_fcol <- function(object,
                          markers,
                          fcol = "markers",
                          type = "l",
                          lty = 1,
                          fractions = sampleNames(object),
                          ylab = "Intensity",
                          xlab = "Fractions",
                          ylim,
                          ...) {
    .data <- exprs(object)
    unkncol <- getUnknowncol()
    stockcol <- getStockcol()
    fData(object)[, fcol] <- factor(fData(object)[, fcol])
    lvs <- levels(fData(object)[, fcol])
    if ("unknown" %in% lvs) {
        i <- which(lvs == "unknown")
        lvs <- c(lvs[-i], lvs[i])
        fData(object)[, fcol] <- factor(fData(object)[, fcol],
                                        levels = lvs)
    }
    ukn <- fData(object)[, fcol] == "unknown"
    .fcol <- fData(object)[, fcol]
    col <- stockcol[as.numeric(.fcol)]
    col[ukn] <- unkncol
    matlines(t(.data),
             col = col,
             type = type,
             lty = lty,
             ...)
}


##' Produces a line plot showing the feature abundances
##' across the fractions.
##'
##' @title Plots the distribution of features across fractions
##' @param object An instance of class \code{MSnSet}.
##' @param markers A \code{character}, \code{numeric} or
##'     \code{logical} of appropriate length and or content used to
##'     subset \code{object} and define the organelle markers.
##' @param fcol Feature meta-data label (fData column name) defining
##'     the groups to be differentiated using different colours. If
##'     \code{NULL} (default) ignored and \code{mcol} and \code{pcol}
##'     are used.
##' @param mcol A \code{character} define the colour of the marker
##'     features.  Default is \code{"steelblue"}.
##' @param pcol A \code{character} define the colour of the
##'     non-markers features.  Default is the colour used for features
##'     of unknown localisation, as returned by
##'     \code{\link{getUnknowncol}}.
##' @param alpha A numeric defining the alpha channel (transparency)
##'     of the points, where \code{0 <= alpha <= 1}, 0 and 1 being
##'     completely transparent and opaque.
##' @param type Character string defining the type of lines. For
##'     example \code{"p"} for points, \code{"l"} for lines,
##'     \code{"b"} for both. See \code{plot} for all possible types.
##' @param lty Vector of line types for the marker profiles. Default
##'     is 1 (solid). See \code{\link{par}} for details.
##' @param fractions A \code{character} defining the \code{phenoData}
##'     variable to be used to label the fraction along the x
##'     axis. Default is to use \code{sampleNames(object)}.
##' @param ylab y-axis label. Default is "Intensity".
##' @param xlab x-axis label. Default is "Fractions".
##' @param ylim A numeric vector of length 2, giving the y coordinates
##'     range.
##' @param ... Additional parameters passed to \code{\link{plot}}.
##' @return Used for its side effect of producing a feature
##'     distribution plot. Invisibly returns the data matrix.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(tan2009r1)
##' j <- which(fData(tan2009r1)$markers == "mitochondrion")
##' i <- which(fData(tan2009r1)$PLSDA == "mitochondrion")
##' plotDist(tan2009r1[i, ], markers = featureNames(tan2009r1)[j])
##' plotDist(tan2009r1[i, ], markers = featureNames(tan2009r1)[j],
##'          fractions = "Fractions")
##' ## plot and colour all marker profiles
##' tanmrk <- markerMSnSet(tan2009r1)
##' plotDist(tanmrk, fcol = "markers")
plotDist <- function(object,
                     markers,
                     fcol = NULL,
                     mcol = "steelblue",
                     pcol = getUnknowncol(),
                     alpha = 0.3,
                     type = "b",
                     lty = 1,
                     fractions = sampleNames(object),
                     ylab = "Intensity",
                     xlab = "Fractions",
                     ylim,
                     ...) {
    .data <- exprs(object)
    if (missing(ylim))
        ylim <- range(.data)
    n <- nrow(.data)
    m <- ncol(.data)
    if (is.character(fractions) & length(fractions) == 1) {
        if (sum(fractions %in% names(pData(object))) != 1)
            stop("'fractions' must be a single pData name.")
        fractions <- as.character(pData(object)[, fractions])
    }
    plot(0, ylim = ylim, xlim = c(1, m),
         xlab = xlab, ylab = ylab,
         type = "n", xaxt = "n")
    axis(1, at = seq_len(m), labels = fractions)
    if (!is.null(fcol)) {
        ## plot and colour all profiles according to fcol
        plotDist_fcol(object, fcol = fcol, ...)
    } else {
        pcol <- col2hcl(pcol, alpha = alpha)
        matlines(t(.data),
                 lty = "solid",
                 col = pcol)
        if (!missing(markers)) {
            mcol <- col2hcl(mcol)
            .mrk <- exprs(object[markers, ])
            matlines(t(.mrk),
                     col = mcol,
                     type = type,
                     lty = lty,
                     ...)
        }
    }
    invisible(t(.data))
}


##' The function checks if the number of clusters to be highlighted is
##' larger than the number of available colours and plotting
##' characters and returns the necessary variables to use both to
##' distinguish all clusters.
##'
##' @title Are there many clusters to plot?
##' @param object An \code{MSnSet} instance
##' @param fcol Feature variable name used as markers
##' @param stockcol Vector of colours.
##' @param stockpch Vector of plotting characters.
##' @return A \code{list} of necessary variables for plot and legend
##' printing. See code for details.
##' @author Laurent Gatto
##' @noRd
.isbig <- function(object, fcol, stockcol, stockpch) {
    if (is.null(fcol))
        return(list(big = FALSE, toobig = FALSE,
                    nclst = NA,
                    ncol = NA, npch = NA,
                    k = NA, kk = NA, jj = NA))
    ## number of clusters to be coloured
    clsts <- unique(fData(object)[, fcol])
    nclst <- length(clsts[clsts != "unknown"])
    ## number of available colours
    ncol <- length(stockcol)
    ## number of available plotting characters
    npch <- length(stockpch)
    big <- (nclst > ncol)
    toobig <- (nclst > ncol * npch)
    if (big) {
        if (toobig) warning(paste0("Not enough colours and pch.\n",
                                   "Some classes will not be coloured."),
                            call. = FALSE)
        else message("Not enough colours: using colours and pch.")
    }
    k <- nclst / ncol
    ## pch indices
    kk <- rep(1:ceiling(k), each = ncol)[1:nclst]
    ## colour indices
    jj <- rep(1:ncol, ceiling(k))[1:nclst]
    return(list(big = big, toobig = toobig,
                nclst = nclst,
                ncol = ncol, npch = npch,
                k = k, kk = kk, jj = jj))
}


## Available pRoloc visualisation methods
pRolocVisMethods <- c("PCA", "MDS", "kpca", "lda", "t-SNE", "nipals",
                      "hexbin", "none")

## Available plot2D methods
plot2Dmethods <- c(pRolocVisMethods, "scree")

##' Plot organelle assignment data and results.
##'
##' Generate 2 or 3 dimensional feature distribution plots to
##' illustrate localistation clusters. Rows/features containing
##' \code{NA} values are removed prior to dimension reduction except
##' for the \code{"nipals"} method. For this method, it is advised to
##' set the method argument `ncomp` to a low number of dimensions to
##' avoid computing all components when analysing large datasets.
##'
##' \code{plot3D} relies on the ##' \code{rgl} package, that will be
##' loaded automatically.
##'
##' \itemize{
##'
##' \item Note that \code{plot2D} has been update in version 1.3.6 to
##'        support more organelle classes than colours defined in
##'        \code{\link{getStockcol}}. In such cases, the default
##'        colours are recycled using the default plotting characters
##'        defined in \code{\link{getStockpch}}. See the example for
##'        an illustration. The \code{alpha} argument is also
##'        depreciated in version 1.3.6. Use \code{setStockcol} to set
##'        colours with transparency instead. See example below.
##'
##' \item Version 1.11.3: to plot data as is, i.e. without any
##'       transformation, \code{method} can be set to "none" (as
##'       opposed to passing pre-computed values to \code{method} as a
##'       \code{matrix}, in previous versions). If \code{object} is an
##'       \code{MSnSet}, the untransformed values in the assay data
##'       will be plotted. If \code{object} is a \code{matrix} with
##'       coordinates, then a matching \code{MSnSet} must be passed to
##'       \code{methargs}.
##' }
##'
##'
##' @param object An instance of class \code{MSnSet}.
##' @param fcol Feature meta-data label (fData column name) defining
##'     the groups to be differentiated using different
##'     colours. Default is \code{markers}. Use \code{NULL} to
##'     suppress any colouring.
##' @param fpch Featre meta-data label (fData column name) desining
##'     the groups to be differentiated using different point symbols.
##' @param unknown A \code{character} (default is \code{"unknown"})
##'     defining how proteins of unknown/un-labelled localisation are
##'     labelled.
##' @param dims A \code{numeric} of length 2 (or 3 for \code{plot3D})
##'     defining the dimensions to be plotted. Defaults are \code{c(1,
##'     2)} and \code{c(1, 2, 3)}.  Always \code{1:2} for MDS.
##' @param score A numeric specifying the minimum organelle assignment
##'     score to consider features to be assigned an organelle. (not
##'     yet implemented).
##' @param method A \code{character} describe how to transform the
##'     data or what to plot. One of \code{"PCA"} (default),
##'     \code{"MDS"}, \code{"kpca"}, \code{"nipals"}, \code{"t-SNE"}
##'     or \code{"lda"}, defining what dimensionality reduction is
##'     applied: principal component analysis (see
##'     \code{\link{prcomp}}), classical multidimensional scaling (see
##'     \code{\link{cmdscale}}), kernel PCA (see
##'     \code{\link[kernlab]{kpca}}), nipals (principal component
##'     analysis by NIPALS, non-linear iterative partial least squares
##'     which support missing values; see
##'     \code{\link[nipals]{nipals}}) t-SNE (see
##'     \code{\link[Rtsne]{Rtsne}}) or linear discriminant analysis
##'     (see \code{\link[MASS]{lda}}). The last method uses
##'     \code{fcol} to defined the sub-cellular clusters so that the
##'     ration between within ad between cluster variance is
##'     maximised. All the other methods are unsupervised and make use
##'     \code{fcol} only to annotate the plot. Prior to t-SNE,
##'     duplicated features are removed and a message informs the user
##'     if such filtering is needed.
##'
##'     \code{"scree"} can also be used to produce a scree
##'     plot. \code{"hexbin"} applies PCA to the data and uses
##'     bivariate binning into hexagonal cells from
##'     \code{\link[hexbin]{hexbin}} to emphasise cluster density.
##'
##'     If none is used, the data is plotted as is, i.e. without any
##'     transformation. In this case, \code{object} can either be an
##'     \code{MSnSet} or a \code{matrix} (as invisibly returned by
##'     \code{plot2D}). This enables to re-generate the figure without
##'     computing the dimensionality reduction over and over again,
##'     which can be time consuming for certain methods. If \code{object}
##'     is a \code{matrix}, an \code{MSnSet} containing the feature
##'     metadata must be provided in \code{methargs} (see below for
##'     details).
##'
##'     Available methods are listed in \code{plot2Dmethods}.
##'
##' @param methargs A \code{list} of arguments to be passed when
##'     \code{method} is called. If missing, the data will be scaled
##'     and centred prior to PCA and t-SNE (i.e. \code{Rtsne}'s
##'     arguments \code{pca_center} and \code{pca_scale} are set to
##'     \code{TRUE}). If \code{method = "none"} and \code{object} is a
##'     \code{matrix}, then the first and only argument of
##'     \code{methargs} must be an \code{MSnSet} with matching
##'     features with \code{object}.
##' @param axsSwitch A \code{logical} indicating whether the axes
##'     should be switched.
##' @param mirrorX A \code{logical} indicating whether the x axis
##'     should be mirrored?
##' @param mirrorY A \code{logical} indicating whether the y axis
##'     should be mirrored?
##' @param col A \code{character} of appropriate length defining
##'     colours.
##' @param pch A \code{character} of appropriate length defining point
##'     character.
##' @param cex Character expansion.
##' @param index A \code{logical} (default is \code{FALSE}, indicating
##'     of the feature indices should be plotted on top of the
##'     symbols.
##' @param idx.cex A \code{numeric} specifying the character expansion
##'     (default is 0.75) for the feature indices. Only relevant when
##'     \code{index} is TRUE.
##' @param addLegend A character indicating where to add the
##'     legend. See \code{\link{addLegend}}for details. If missing
##'     (default), no legend is added.
##' @param identify A logical (default is \code{TRUE}) defining if
##'     user interaction will be expected to identify individual data
##'     points on the plot. See also \code{\link{identify}}.
##' @param plot A \code{logical} defining if the figure should be
##'     plotted.  Useful when retrieving data only. Default is
##'     \code{TRUE}.
##' @param grid A \code{logical} indicating whether a grid should
##'     be plotted. Default is \code{TRUE}.
##' @param ... Additional parameters passed to \code{plot} and
##'     \code{points}.
##' @return Used for its side effects of generating a plot.  Invisibly
##'     returns the 2 or 3 dimensions that are plotted.
##' @author Laurent Gatto <lg390@@cam.ac.uk>
##' @seealso \code{\link{addLegend}} to add a legend to \code{plot2D}
##'     figures (the legend is added by default on \code{plot3D}) and
##'     \code{\link{plotDist}} for alternative graphical
##'     representation of quantitative organelle proteomics
##'     data. \code{\link{plot2Ds}} to overlay 2 data sets on the same
##'     PCA plot. The \code{\link{plotEllipse}} function can be used
##'     to visualise TAGM models on PCA plots with ellipses.
##' @aliases plot2Dmethods
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' plot2D(dunkley2006, fcol = NULL)
##' plot2D(dunkley2006, fcol = NULL, col = "black")
##' plot2D(dunkley2006, fcol = "markers")
##' addLegend(dunkley2006,
##'           fcol = "markers",
##'           where = "topright",
##'           cex = 0.5, bty = "n", ncol = 3)
##' title(main = "plot2D example")
##' ## available methods
##' plot2Dmethods
##' plot2D(dunkley2006, fcol = NULL, method = "kpca", col = "black")
##' plot2D(dunkley2006, fcol = NULL, method = "kpca", col = "black",
##'        methargs = list(kpar = list(sigma = 1)))
##' plot2D(dunkley2006, method = "lda")
##' plot2D(dunkley2006, method = "hexbin")
##' ## Using transparent colours
##' setStockcol(paste0(getStockcol(), "80"))
##' plot2D(dunkley2006, fcol = "markers")
##' ## New behavious in 1.3.6 when not enough colours
##' setStockcol(c("blue", "red", "green"))
##' getStockcol() ## only 3 colours to be recycled
##' getMarkers(dunkley2006)
##' plot2D(dunkley2006)
##' ## reset colours
##' setStockcol(NULL)
##' plot2D(dunkley2006, method = "none") ## plotting along 2 first fractions
##' plot2D(dunkley2006, dims = c(3, 5), method = "none") ## plotting along fractions 3 and 5
##' ## pre-calculate PC1 and PC2 coordinates
##' pca <- plot2D(dunkley2006, plot=FALSE)
##' head(pca)
##' plot2D(pca, method = "none", methargs  = list(dunkley2006))
plot2D <- function(object,
                   fcol = "markers",
                   fpch,
                   unknown = "unknown",
                   dims = 1:2,
                   score = 1, ## TODO
                   method = "PCA",
                   methargs,
                   axsSwitch = FALSE,
                   mirrorX = FALSE,
                   mirrorY = FALSE,
                   col,
                   pch,
                   cex,
                   index = FALSE,
                   idx.cex = 0.75,
                   addLegend,
                   identify = FALSE,
                   plot = TRUE,
                   grid = TRUE,
                   ...) {
    ## handling deprecated outliers argument
    a <- as.list(match.call()[-1])
    if ("outliers" %in% names(a))
        stop("'outliers' is deprecated. Use xlim/ylim to focus your plot")
    if (!missing(col)) {
        stockcol <- col
    } else {
        stockcol <- getStockcol()
    }
    if (!missing(pch)) {
        stockpch <- pch
        userpch <- TRUE
    } else {
        stockpch <- getStockpch()
        userpch <- FALSE
    }
    unknowncol <- getUnknowncol()
    unknownpch <- getUnknownpch()
    method <- match.arg(method, plot2Dmethods)

    if (length(dims) > 2) {
        warning("Using first two dimensions of ", dims)
        dims <- dims[1:2]
    }
    k <- max(dims)
    if (any(is.na(object)) & method != "nipals") {
        n0 <- nrow(object)
        object <- filterNA(object)
        n1 <- nrow(object)
        if (n1 == 0)
            stop("No rows left after removing NAs!")
        else
            message("Removed ", n0 - n1, " row(s) with 'NA' values.\n",
                    "Consider using 'nipals' to retain all features.")
    }
    if (method == "hexbin") {
        requireNamespace("hexbin")
        if (missing(methargs))
            methargs <- list(scale = TRUE, center = TRUE)
        .data <- plot2D(object,
                        method = "PCA",
                        fcol = NULL,
                        dims = dims,
                        methargs = methargs,
                        plot = FALSE)
        cramp <- colorRampPalette(c("grey90", "blue"))
        plot(hexbin::hexbin(.data), colramp = cramp, ...)
        return(invisible(.data))
    } else if (method == "scree") {
        if (missing(methargs))
            methargs <- list(scale = TRUE, center = TRUE)
        .pca <- do.call(prcomp, c(list(x = exprs(object)),
                                  methargs))
        ## plot(.pca, npcs = ncol(.data))
        .vars <- (.pca$sdev)^2
        if (plot)
            barplot(.vars / sum(.vars) * 100,
                    ylab = "Percentage of total variance",
                    xlab = paste0("PC 1 to ", ncol(object)))
        .data <- .vars / sum(.vars) * 100
        plot <- FALSE
        fcol <- NULL
    } else if (method == "lda") {
        if (!is.null(fcol) && !fcol %in% fvarLabels(object))
            stop("'", fcol, "' not found in feature variables.")
        requireNamespace("MASS")
        X <- data.frame(exprs(object))
        gr <- getMarkers(object, fcol = fcol, verbose = FALSE)
        train <- which(gr != "unknown")
        if (missing(methargs)) {
            z <- MASS::lda(X, grouping = gr, subset = train)
        } else {
            z <- do.call(MASS::lda, c(list(x = X,
                                           grouping = gr,
                                           subset = train),
                                      methargs))
        }
        p <- predict(z, X)
        .data <- p$x[, dims]
        tr <- round(z$svd^2 / sum(z$svd^2), 4L) * 100
        .xlab <- paste0("LD", dims[1], " (", tr[dims[1]], "%)")
        .ylab <- paste0("LD", dims[2], " (", tr[dims[2]], "%)")
    } else if (method == "t-SNE") {
        if (!requireNamespace("Rtsne") && packageVersion("Rtsne") >= 0.13)
            stop("Please install the Rtsne (>= 0.13) package to make use if this functionality.")
        if (missing(methargs))
            methargs <- list(pca_scale = TRUE, pca_center = TRUE)
        e <- exprs(object)
        nr0 <- nrow(e)
        e <- unique(e)
        if (nrow(e) < nr0) {
            message("Only keeping unique features, dropped ", nr0 - nrow(e), ".")
            object <- object[rownames(e), ] ## see #108
        }
        .data <- do.call(Rtsne::Rtsne,
                         c(list(X = e),
                           dims = k,
                           methargs))
        .data <- .data$Y[, dims]
        .xlab <- paste("Dimension", dims[1])
        .ylab <- paste("Dimension", dims[2])
        colnames(.data) <- c(.xlab, .ylab)
        rownames(.data) <- rownames(e)
    } else if (method == "PCA") {
        if (missing(methargs))
            methargs <- list(scale = TRUE, center = TRUE)
        .data <- dimred(object, method = method, methargs = methargs)
        .data <- .data[, dims]
        .xlab <- colnames(.data)[1]
        .ylab <- colnames(.data)[2]
    } else if (method == "nipals") {
        if (!requireNamespace("nipals"))
            stop("Please install the nipals package to make use if this functionality.")
        if (missing(methargs))
            methargs <- list(scale = TRUE, center = TRUE, ncomp = ncol(object))
        .data <- dimred(object, method = method, methargs = methargs)
        .data <- .data[, dims]
        .xlab <- colnames(.data)[1]
        .ylab <- colnames(.data)[2]
    } else if (method == "MDS")  {
        if (!missing(methargs))
            warning("'methargs' ignored for MDS")
        ## TODO - use other distances
        .data <- cmdscale(dist(exprs(object),
                               method = "euclidean",
                               diag = FALSE,
                               upper = FALSE),
                          k = 2)
        .xlab <- paste("Dimension 1")
        .ylab <- paste("Dimension 2")
        colnames(.data) <- c(.xlab, .ylab)
    } else if (method == "kpca") {
        if (missing(methargs)) {
            .kpca <- kpca(exprs(object))
        } else {
            .kpca <- do.call(kpca, c(list(x = exprs(object)),
                                     methargs))
        }
        .data <- rotated(.kpca)[, dims]
        .vars <- (eig(.kpca)/sum(eig(.kpca)))[dims]
        .vars <- round(100 * .vars, 2)
        .xlab <- paste0("PC", dims[1], " (", .vars[1], "%)")
        .ylab <- paste0("PC", dims[2], " (", .vars[2], "%)")
        colnames(.data) <- c(.xlab, .ylab)
    } else { ## none
        if (inherits(object, "MSnSet")) {
            .data <- exprs(object)[, dims]
        } else if (is.matrix(object)) {
            .data <- object[, dims]
            ## we still need the feature variables
            object <- methargs[[1]]
            if (!inherits(object, "MSnSet"))
                stop(paste("If method == \"none\", and object is a 'matrix',",
                           "the feature metadata must be provided as an 'MSnSet'",
                           "(the object matching the coordinate matrix) in 'methargs'"))
            if (nrow(.data) != nrow(object))
                stop("Number of features in the matrix and feature metadata differ.")
            if (!all.equal(rownames(.data), featureNames(object)))
                warning("Matrix rownames and feature names don't match")
        } else stop("object must be an 'MSnSet' or a 'matrix' (if method == \"none\").")
        .xlab <- colnames(.data)[1]
        .ylab <- colnames(.data)[2]
    }
    ## Now, object must be an MSnSet - if it was a matrix, it has been
    ## replaced by methargs[[1]] (see method = "none" above).
    stopifnot(inherits(object, "MSnSet"))
    if (isMrkMat(object, fcol))
        stop("To visualise a marker matrix, use 'pRolocVis' from package pRolocGUI (>= 1.5.2).")
    if (!is.null(fcol)) stopifnot(isMrkVec(object, fcol))
    if (!is.null(fcol) && !fcol %in% fvarLabels(object))
        stop("'", fcol, "' not found in feature variables.")
    if (!missing(fpch) && !fpch %in% fvarLabels(object))
        stop("'", fpch, "' not found in feature variables.")
    ## mirror irrespective of plotting
    if (mirrorX)
        .data[, 1] <- -.data[, 1]
    if (mirrorY)
        .data[, 2] <- -.data[, 2]
    if (plot) {
        if (axsSwitch) {
            .data <- .data[, 2:1]
            .tmp <- .xlab
            .xlab <- .ylab
            .ylab <- .tmp
        }
        if (missing(col))
            col <- rep(unknowncol, nrow(.data))
        if (missing(pch))
            pch <- rep(unknownpch, nrow(.data))
        if (missing(cex)) {
            cex <- rep(1, nrow(.data))
        } else {
            if (length(cex) == 1) {
                cex <- rep(cex, nrow(.data))
            } else {
                .n <- nrow(.data) %/% length(cex)
                .m <- nrow(.data) %% length(cex)
                cex <- c(rep(cex, .n),
                         cex[.m])
            }
        }
        stopifnot(length(cex) == nrow(.data))
        if (!is.null(fcol)) {
            nullfcol <- FALSE
            fData(object)[, fcol] <- factor(fData(object)[, fcol])
            lvs <- levels(fData(object)[, fcol])
            if ("unknown" %in% lvs) {
                i <- which(lvs == "unknown")
                lvs <- c(lvs[-i], lvs[i])
                fData(object)[, fcol] <- factor(fData(object)[, fcol],
                                                levels = lvs)
            }
            ukn <- fData(object)[, fcol] == unknown
            .fcol <- fData(object)[, fcol]
            col <- stockcol[as.numeric(.fcol)]
            col[ukn] <- unknowncol
        } else {
            nullfcol <- TRUE
            ukn <- rep(TRUE, nrow(.data))
        }
        if (!missing(fpch)) {
            .fpch <- factor(fData(object)[, fpch])
            pch <- stockpch[as.numeric(.fpch)]
        } else {
            pch <- rep(stockpch[1], nrow(.data))
        }
        ## don't set this if fcol is null and pch was passed, as then
        ## we want all points to be plotted with the user-defined pch
        if (!(nullfcol & userpch))
            pch[ukn] <- unknownpch
        isbig <- .isbig(object, fcol, stockcol, stockpch)

        if (is.null(fcol)) {
            plot(.data, xlab = .xlab, ylab = .ylab, col = col,
                 pch = pch, cex = cex, ...)
        } else if (isbig[["big"]]) {
            plot(.data, xlab = .xlab, ylab = .ylab,
                 type = "n", ...)
            points(.data[ukn, 1], .data[ukn, 2],
                   col = col[ukn],
                   pch = pch[ukn], cex = cex[ukn], ...)
            clst <- levels(factor(fData(object)[, fcol]))
            clst <- clst[clst != "unknown"]
            for (i in 1:isbig[["nclst"]]) {
                sel <- fData(object)[, fcol] == clst[i]
                points(.data[sel, 1], .data[sel, 2],
                       cex = cex[!ukn],
                       col = stockcol[isbig[["jj"]][i]],
                       pch = stockpch[isbig[["kk"]][i]])
            }
        } else {
            plot(.data, xlab = .xlab, ylab = .ylab,
                 type = "n", ...)
            points(.data[ukn, 1], .data[ukn, 2],
                   col = col[ukn],
                   pch = pch[ukn], cex = cex[ukn])
            points(.data[!ukn, 1], .data[!ukn, 2],
                   col = col[!ukn],
                   pch = pch[!ukn], cex = cex[!ukn])
        }
        if (!missing(addLegend)) {
            where <- match.arg(addLegend,
                               c("bottomright", "bottom",
                                 "bottomleft", "left", "topleft",
                                 "top", "topright", "right", "center"))
            addLegend(object, fcol = fcol, where = where)
        }
        if (grid) grid()
        if (index) {
            text(.data[, 1], .data[, 2], 1:nrow(.data), cex = idx.cex)
        }
        if (identify) {
            ids <- identify(.data[, 1], .data[, 2],
                            rownames(.data))
            return(ids)
        }
    }
    invisible(.data)
}


##' Adds a legend to a \code{\link{plot2D}} figure.
##'
##' The function has been updated in version 1.3.6 to recycle the
##' default colours when more organelle classes are provided. See
##' \code{\link{plot2D}} for details.
##'
##' @title Adds a legend
##' @param object An instance of class \code{MSnSet}
##' @param fcol Feature meta-data label (fData column name) defining
##'     the groups to be differentiated using different
##'     colours. Default is \code{markers}.
##' @param where One of \code{"bottomleft"} (default),
##'     \code{"bottomright"}, \code{"topleft"}, \code{"topright"} or
##'     \code{"other"} defining the location of the
##'     legend. \code{"other"} opens a new graphics device, while the
##'     other locations are passed to \code{\link{legend}}.
##' @param col A \code{character} defining point colours.
##' @param bty Box type, as in \code{legend}. Default is set to
##'     \code{"n"}.
##' @param ... Additional parameters passed to \code{\link{legend}}.
##' @return Invisibly returns \code{NULL}
##' @author Laurent Gatto
addLegend <- function(object,
                      fcol = "markers",
                      where = c("bottomleft", "bottom", "bottomright",
                                "left", "topleft", "top", "topright",
                                "right", "center", "other"),
                      col,
                      bty = "n",
                      ...) {
    where <- match.arg(where)
    fData(object)[, fcol] <- as.factor(fData(object)[, fcol])
    lvs <- levels(fData(object)[, fcol])
    if ("unknown" %in% lvs) {
        i <- which(lvs == "unknown")
        lvs <- c(lvs[-i], lvs[i])
        fData(object)[, fcol] <- factor(fData(object)[, fcol],
                                        levels = lvs)
    } else {
        fData(object)[, fcol] <- factor(fData(object)[, fcol])
    }

    if (where == "other") {
        dev.new()
        where <- "center"
        plot(0, type = "n", bty = "n",
             xaxt = "n", yaxt = "n",
             xlab = "", ylab = "")
    }
    if (is.null(fcol))
        fcol <- "markers"
    if (!fcol %in% fvarLabels(object))
        stop("'", fcol, "' not found in feature variables.")
    if (missing(col)) {
        stockcol <- getStockcol()
    } else {
        stockcol <- col
    }
    stockpch <- getStockpch()
    unknowncol <- getUnknowncol()
    unknownpch <- getUnknownpch()
    txt <- levels(as.factor(fData(object)[, fcol]))
    isbig <- .isbig(object, fcol, stockcol, stockpch)
    if (isbig[["big"]]) {
        col <- stockcol[isbig[["jj"]]]
        pch <- stockpch[isbig[["kk"]]]
        if ("unknown" %in% txt) {
            col <- c(col, unknowncol)
            pch <- c(pch, unknownpch)
        }
    } else {
        col <- stockcol[seq_len(length(txt))]
        pch <- rep(stockpch[1], length(txt))
        if ("unknown" %in% txt) {
            i <- which(txt == "unknown")
            col[-i] <- col[-length(col)]
            col[i] <- unknowncol
            pch[-i] <- pch[-length(pch)]
            pch[i] <- unknownpch
        }
    }
    if ("bty" %in% names(pairlist(...))) legend(where, txt, col = col, pch = pch, ...)
    else legend(where, txt, col = col, pch = pch, bty = bty, ...)
    invisible(NULL)
}

##' Highlights a set of features of interest given as a
##' \code{FeaturesOfInterest} instance on a PCA plot produced by
##' code{plot2D} or \code{plot3D}. If none of the features of interest
##' are found in the \code{MSnset}'s \code{featureNames}, an warning
##' is thrown.
##'
##' @title Highlight features of interest on a spatial proteomics plot
##'
##' @param object The main dataset described as an \code{MSnSet} or a
##'     \code{matrix} with the coordinates of the features on the PCA
##'     plot produced (and invisibly returned) by \code{plot2D}.
##'
##' @param foi An instance of \code{\linkS4class{FeaturesOfInterest}},
##'     or, alternatively, a \code{character} of feautre names.
##'
##' @param labels A \code{character} of length 1 with a feature
##'     variable name to be used to label the features of
##'     interest. This is only valid if \code{object} is an
##'     \code{MSnSet}. Alternatively, if \code{TRUE}, then
##'     \code{featureNames(object)} (or code{rownames(object)}, if
##'     \code{object} is a \code{matrix}) are used. Default is
##'     missing, which does not add any label.s
##'
##' @param args A named list of arguments to be passed to
##'     \code{plot2D} if the PCA coordinates are to be
##'     calculated. Ignored if the PCA coordinates are passed
##'     directly, i.e. \code{object} is a \code{matrix}.
##'
##' @param ... Additional parameters passed to \code{points} or
##'     \code{text} (when \code{labels} is \code{TRUE}) when adding to
##'     \code{plot2D}, or \code{spheres3d} or \code{text3d} when
##'     adding the \code{plot3D}
##' @return NULL; used for its side effects.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data("tan2009r1")
##' x <- FeaturesOfInterest(description = "A test set of features of interest",
##'                         fnames = featureNames(tan2009r1)[1:10],
##'                         object = tan2009r1)
##'
##' ## using FeaturesOfInterest or feature names
##' par(mfrow = c(2, 1))
##' plot2D(tan2009r1)
##' highlightOnPlot(tan2009r1, x)
##' plot2D(tan2009r1)
##' highlightOnPlot(tan2009r1, featureNames(tan2009r1)[1:10])
##'
##' .pca <- plot2D(tan2009r1)
##' head(.pca)
##' highlightOnPlot(.pca, x, col = "red")
##' highlightOnPlot(tan2009r1, x, col = "red", cex = 1.5)
##' highlightOnPlot(tan2009r1, x, labels = TRUE)
##'
##' .pca <- plot2D(tan2009r1, dims = c(1, 3))
##' highlightOnPlot(.pca, x, pch = "+", dims = c(1, 3))
##' highlightOnPlot(tan2009r1, x, args = list(dims = c(1, 3)))
##'
##' .pca2 <- plot2D(tan2009r1, mirrorX = TRUE, dims = c(1, 3))
##' ## previous pca matrix, need to mirror X axis
##' highlightOnPlot(.pca, x, pch = "+", args = list(mirrorX = TRUE))
##' ## new pca matrix, with X mirrors (and 1st and 3rd PCs)
##' highlightOnPlot(.pca2, x, col = "red")
##'
##' plot2D(tan2009r1)
##' highlightOnPlot(tan2009r1, x)
##' highlightOnPlot(tan2009r1, x, labels = TRUE, pos = 3)
##' highlightOnPlot(tan2009r1, x, labels = "Flybase.Symbol", pos = 1)
highlightOnPlot <- function(object, foi, labels, args = list(), ...) {
    if (is.character(foi))
        foi <- FeaturesOfInterest(description = "internally created",
                                  fnames = foi)
    stopifnot(inherits(foi, "FeaturesOfInterest"))

    if (!fnamesIn(foi, object)) {
        warning("None of the features of interest are present in the data.")
        return(invisible(NULL))
    }
    if (inherits(object, "MSnSet")) {
        .args <- list(object = object, plot = FALSE, fcol = NULL)
        args <- c(args, .args)
        .pca <- do.call(plot2D, args = args)
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
        text(.pca[sel, 1], .pca[sel, 2], labels[sel], ...)
    } else {
        points(.pca[sel, 1], .pca[sel, 2], ...)
    }
}


## Tests whether the object is visualisation method available
## in pRoloc
.validpRolocVisMethod <- function(object) {
    if (class(object) == "matrix" && ncol(object) == 2)
        return(TRUE)
    else
        if (object %in% pRolocVisMethods)
            return(TRUE)
    else
        return(FALSE)
}

##' The function calculates and draws convex hulls on top of a plot
##' produced by \code{\link{plot2D}}. The hulls are computed in the
##' two componnets as defined by the \code{plot2D} method and
##' dimensions.
##'
##' @title Overlay convex hulls onto a \code{plot2D} plot
##' @param object An instance of class \code{MSnSet}
##' @param fcol A \code{character} that defined which feature data
##'     variable to be used as marker definition. Default is
##'     \code{"markers"}.
##' @param ... Additional parameters passed to `plot2D`.
##' @return Invisibly returns \code{NULL}.
##' @author Laurent Gatto
##' @noRd
##' @examples
##' library("pRolocdata")
##' data(E14TG2aS1)
##' plot2D(E14TG2aS1)
##' addConvexHulls(E14TG2aS1)
addConvexHulls <- function(object, fcol = "markers", ...) {
    X <- plot2D(object, fcol = fcol, ...)
    mm <- getMarkerClasses(object, fcol = fcol)
    for (.mm in mm) {
        sel <- getMarkers(object, fcol = fcol, verbose = FALSE) == .mm
        .X <- X[sel, ]
        hpts <- chull(.X)
        hpts <- c(hpts, hpts[1])
        lines(.X[hpts, ])
    }
    invisible(NULL)
}


getMarkerCols <- function(object, fcol = "markers") {
    stopifnot(fcol %in% fvarLabels(object))
    .fcol <- factor(fData(object)[, fcol])
    lvs <- levels(.fcol)
    if ("unknown" %in% lvs) {
        i <- which(lvs == "unknown")
        lvs <- c(lvs[-i], lvs[i])
        .fcol <- factor(.fcol, levels = lvs)
    }
    ukn <- .fcol == "unknown"
    col <- getStockcol()[as.numeric(.fcol)]
    col[ukn] <- getUnknowncol()
    col
}


##' The function plots marker consensus profiles obtained from mrkConsProfile
##'
##' @title Plot marker consenses profiles.
##' @param object A `matrix` containing marker consensus profiles as output from
##'     [mrkConsProfiles()].
##' @param order Order for markers (optional).
##' @param plot A `logical(1)` defining whether the heatmap should be plotted.
##'     Default is `TRUE`.
##' @return Invisibly returns `ggplot2` object.
##' @author Tom Smith
##' @md
##' @examples
##' library("pRolocdata")
##' data(E14TG2aS1)
##' hc <- mrkHClust(E14TG2aS1, plot = FALSE)
##' mm <- getMarkerClasses(E14TG2aS1)
##' ord <- levels(factor(mm))[order.dendrogram(hc)]
##' fmat <- mrkConsProfiles(E14TG2aS1)
##' plotConsProfiles(fmat, order = ord)
plotConsProfiles <- function(object, order = NULL, plot = TRUE) {
    feature <- NULL
    fmatlong <- cbind(expand.grid("feature" = rownames(object),
                                  "sample" = colnames(object),
                                  stringsAsFactors = FALSE),
                      "intensity" = as.vector(object))
    if (!is.null(order))
        fmatlong$feature <- factor(fmatlong$feature, order)

    fmatlong$sample <- factor(fmatlong$sample, colnames(object))

    p <- ggplot(fmatlong, aes(sample, feature, fill = intensity)) +
        geom_tile() +
        scale_fill_continuous(low = "white", high = "#56B4E9",
                              limits = c(0, NA),
                              name = "Intensity") +
        theme_bw() +
        xlab("Sample") +
        ylab("") +
        theme(panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              aspect.ratio = 1) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0))

    if (plot)
        print(p)

    invisible(p)
}
