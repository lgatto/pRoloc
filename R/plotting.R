##' Produces a line plot showing the feature abundances
##' across the fractions.
##'
##' @title Plots the distribution of features across fractions
##' @param object An instance of class \code{MSnSet}.
##' @param markers A \code{character}, \code{numeric} or
##'     \code{logical} of appropriate length and or content used to
##'     subset \code{object} and define the organelle markers.
##' @param mcol A \code{character} define the colour of the marker
##'     features.  Default is \code{"steelblue"}.
##' @param pcol A \code{character} define the colour of the
##'     non-markers features.  Default is \code{"grey90"}.
##' @param alpha A numeric defining the alpha channel (transparency)
##'     of the points, where \code{0 <= alpha <= 1}, 0 and 1 being
##'     completely transparent and opaque.
##' @param type Character string defining the type of lines. For
##'     example \code{"p"} for points, \code{"l"} for lines,
##'     \code{"b"} for both. See \code{plot} for all possible types.
##' @param lty Vector of line types for the marker profiles. Default
##'     is 1 (solid). See \code{\link{par}} for details.
##' @param fractions An optional \code{character} defining the
##'     \code{phenoData} variable to be used to label the fraction
##'     along the x axis. If missing, the \code{phenoData} variables
##'     are searched for a match to \code{fraction}.  If no match is
##'     found, the fractions are labelled as numericals.
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
##' plotDist(tan2009r1[i, ],
##'          markers = featureNames(tan2009r1)[j])
##' title("Mitochondrion")
plotDist <- function(object,
                     markers,
                     mcol = "steelblue",                     
                     pcol = "grey90",
                     alpha = 0.3,
                     type = "b",
                     lty = 1,
                     fractions,
                     ylim,
                     ...) {
  .data <- exprs(object)
  if (missing(ylim))
      ylim <- range(.data)
  n <- nrow(.data)
  m <- ncol(.data)  
  if (missing(fractions)) {
    .frac <- grep("fraction",
                  casefold(names(pData(object))))
    if (length(.frac) == 1) {
      .frac <- as.character(pData(object)[, .frac])
    } else {
      .frac <- seq_len(m)
    }
  }
  plot(0, ylim = ylim,
       xlim = c(1, m),
       ylab = "Intensity",
       xlab = "Fractions",
       type = "n", xaxt = "n")
  axis(1, at = seq_len(m), labels = .frac)
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
pRolocVisMethods <- c("PCA", "MDS", "kpca", "t-SNE", "none")

## Available plot2D methods
plot2Dmethods <- c(pRolocVisMethods, "scree")

##' Plot organelle assignment data and results.
##' 
##' Generate 2 dimensional or feature distribution plots to illustrate
##' localistation clusters. In \code{plot2D}, rows containing
##' \code{NA} values are removed prior to dimention reduction.
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
##' 
##' @param object An instance of class \code{MSnSet}.
##' @param fcol Feature meta-data label (fData column name) defining
##' the groups to be differentiated using different colours. Default
##' is \code{markers}. Use \code{NULL} to suppress any colouring.
##' @param fpch Featre meta-data label (fData column name) desining
##' the groups to be differentiated using different point symbols.
##' @param unknown A \code{character} (default is \code{"unknown"})
##' defining how proteins of unknown/un-labelled localisation are
##' labelled.
##' @param dims A \code{numeric} of length 2 defining the dimensions
##' to be plotted. Always 1:2 for MDS.
##' @param score A numeric specifying the minimum organelle assignment score
##' to consider features to be assigned an organelle. (not yet implemented).
##' 
##' @param method A \code{character} describe how to transform the
##'     data or what to plot. One of \code{"PCA"} (default),
##'     \code{"MDS"}, \code{"kpca"} or \code{"t-SNE"}, defines what
##'     dimensionality reduction is applied: principal component
##'     analysis (see \code{\link{prcomp}}), classical
##'     multidimensional scaling (see \code{\link{cmdscale}}), kernel
##'     PCA (see \code{kernlab::kpca}) or t-SNE (see
##'     \code{tsne::tsne}). \code{"scree"} can also be used to produce
##'     a scree plot. If none is used, the data is plotted as is,
##'     i.e. without any transformation. In this case, \code{object}
##'     can either be an \code{MSnSet} or a \code{matrix} (as
##'     invisibly returned by \code{plot2D}). This enables to
##'     re-generate the figure without computing the dimensionality
##'     reduction over and over again, which can be time consuming for
##'     certain methods. Available methods are listed in
##'     \code{plot2Dmethods}. If \code{object} is a \code{matrix}, an
##'     \code{MSnSet} containing the feature metadata must be provided
##'     in \code{methargs} (see below for details).
##' 
##' @param methargs A \code{list} of arguments to be passed when
##'     \code{method} is called. If missing, the data will be scaled
##'     and centred prior to PCA. If \code{method = "none"} and
##'     \code{object} is a \code{matrix}, then the first and only
##'     argument of \code{methargs} must be an \code{MSnSet} with
##'     matching features with \code{object}.
##' 
##' @param axsSwitch A \code{logical} indicating whether the axes
##'     should be switched.
##' 
##' @param mirrorX A \code{logical} indicating whether the x axis
##'     should be mirrored?
##' 
##' @param mirrorY A \code{logical} indicating whether the y axis
##'     should be mirrored?
##' 
##' @param col A \code{character} of appropriate length defining colours.
##' @param pch A \code{character} of appropriate length defining point character.
##' @param cex Character expansion.
##' @param index A \code{logical} (default is \code{FALSE}, indicating of the
##' feature indices should be plotted on top of the symbols.
##' @param idx.cex A \code{numeric} specifying the character expansion
##' (default is 0.75) for the feature indices. Only relevant when \code{index}
##' is TRUE.
##' @param identify A logical (default is \code{TRUE}) defining if
##' user interaction will be expected to identify individual data
##' points on the plot. See also \code{\link{identify}}.
##' @param plot A \code{logical} defining if the figure should be plotted.
##' Useful when retrieving data only. Default is \code{TRUE}.
##' @param ... Additional parameters passed to \code{plot} and
##' \code{points}.
##' @return Used for its side effects of generating a plot.
##' Invisibly returns the 2d data.
##' @author Laurent Gatto <lg390@@cam.ac.uk>
##' @seealso \code{\link{addLegend}} to add a legend to \code{plot2D}
##'     figures and \code{\link{plotDist}} for alternative graphical
##'     representation of quantitative organelle proteomics
##'     data. \code{\link{plot2Ds}} to overlay 2 data sets on the same
##'     PCA plot.
##' @aliases plot2Dmethods
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' plot2D(dunkley2006, fcol = NULL)
##' ## available methods
##' plot2Dmethods
##' plot2D(dunkley2006, fcol = NULL, method = "kpca")
##' plot2D(dunkley2006, fcol = NULL, method = "kpca",
##'        methargs = list(kpar = list(sigma = 1)))
##' plot2D(dunkley2006, fcol = "markers")
##' addLegend(dunkley2006,
##'           fcol = "markers",
##'           where = "topright",
##'           cex = 0.5, bty = "n", ncol = 3)
##' title(main = "plot2D example")
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
                   identify = FALSE,
                   plot = TRUE,
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
    } else {
        stockpch <- getStockpch()
    }
    unknowncol <- getUnknowncol()
    unknownpch <- getUnknownpch()
    method <- match.arg(method, plot2Dmethods)

    if (length(dims) > 2) {
        warning("Using first two dimensions of ", dims)
        dims <- dims[1:2]
    }
    k <- max(dims)
    if (any(is.na(object))) {
        n0 <- nrow(object)
        object <- filterNA(object)
        n1 <- nrow(object)
        if (n1 == 0)
            stop("No rows left after removing NAs!")
        else
            warning("Removed ", n0 - n1, " row(s) with 'NA' values.")    
    }
    if (method == "scree") {
        if (missing(methargs))
            methargs <- list(scale = TRUE, center = TRUE)
        .pca <- do.call(prcomp, c(list(x = exprs(object)), 
                                  methargs))
        .data <- .pca$x
        plot(.pca, npcs = ncol(.data))
        plot <- FALSE
    } else if (method == "t-SNE") {
        requireNamespace("tsne")
        if (missing(methargs))
            .data <- tsne::tsne(exprs(object), k = k)
        else .data <- do.call(tsne::tsne,
                              c(list(X = exprs(object)),
                                k = k, 
                                methargs))
        .data <- .data[, dims]
        .xlab <- paste("Dimension 1")
        .ylab <- paste("Dimension 2")
        colnames(.data) <- c(.xlab, .ylab)
        rownames(.data) <- featureNames(object)
    } else if (method == "PCA") {
        if (missing(methargs))
            methargs <- list(scale = TRUE, center = TRUE)
        .pca <- do.call(prcomp, c(list(x = exprs(object)), 
                                  methargs))
        .data <- .pca$x[, dims]
        .vars <- (.pca$sdev)^2
        .vars <- (.vars / sum(.vars))[dims]
        .vars <- round(100 * .vars, 2)
        .xlab <- paste0("PC", dims[1], " (", .vars[1], "%)")
        .ylab <- paste0("PC", dims[2], " (", .vars[2], "%)")
        colnames(.data) <- c(.xlab, .ylab)
    } else if (method == "MDS")  { ## MDS
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
        
        col <- rep(unknowncol, nrow(.data))
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
            ukn <- rep(TRUE, nrow(.data))
        }

        if (!missing(fpch)) {   
            .fpch <- factor(fData(object)[, fpch])
            pch <- stockpch[as.numeric(.fpch)]
        } else {
            pch <- rep(stockpch[1], nrow(.data))
        }
        pch[ukn] <- unknownpch
        isbig <- .isbig(object, fcol, stockcol, stockpch)
        
        if (is.null(fcol)) {
            plot(.data, xlab = .xlab, ylab = .ylab, ...)
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
        grid()
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
##' the groups to be differentiated using different colours. Default
##' is \code{markers}. 
##' @param where One of \code{"other"}, \code{"bottomleft"},
##' \code{"bottomright"}, \code{"topleft"} or \code{"topright"} defining
##' the location of the legend. \code{"other"} opens a new graphics device,
##' while the other locations are passed to \code{\link{legend}}.
##' @param col A \code{character} defining point colours.
##' @param bty Box type, as in \code{legend}. Default is set to
##' \code{"n"}.
##' @param ... Additional parameters passed to \code{\link{legend}}.
##' @return Invisibly returns \code{NULL}
##' @author Laurent Gatto
addLegend <- function(object,
                      fcol = "markers",
                      where = c("other", "bottomright", "bottom",
                                "bottomleft", "left", "topleft",
                                "top", "topright", "right", "center"),
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
##' code{plot2D}. If none of the features of interest are found in the
##' \code{MSnset}'s \code{featureNames}, an warning is thrown.
##'
##' @title Highlight features of interest on a plot2D figure
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
##' @param ... Additional parameters passed to \code{points}.
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
        .args <- list(object = object, plot = FALSE)
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
