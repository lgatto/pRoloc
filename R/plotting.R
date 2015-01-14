##' Plot organelle assignment data and results.
##'
##' This is the documentation for the pre-v 1.3.6 function.
##' 
##' Generate 2 dimensional or feature distribution plots to illustrate
##' localistation clusters. In \code{plot2D_v1}, rows containing
##' \code{NA} values are removed prior to dimention reduction.
##' 
##' @param object An instance of class \code{MSnSet}.
##' @param fcol Feature meta-data label (fData column name) defining
##' the groups to be differentiated using different colours. Default
##' is \code{markers}. Use \code{NULL} to suppress any colouring.
##' @param fpch Featre meta-data label (fData column name) desining
##' the groups to be differentiated using different point symbols.
##' @param unknown A \code{character} (default is \code{"unknown"})
##' defining how proteins of unknown localisation are labelled.
##' @param dims A \code{numeric} of length 2 defining the dimensions
##' to be plotted. Always 1:2 for MDS.
##' @param alpha A numeric defining the alpha channel (transparency)
##' of the points, where \code{0 <= alpha <= 1}, 0 and 1 being completely
##' transparent and opaque.
##' @param score A numeric specifying the minimum organelle assignment score
##' to consider features to be assigned an organelle. (not yet implemented).
##' @param outliers A logical specifying whether outliers should be plotted
##' or ignored (default is TRUE, i.e. all points are plotted). Useful when 
##' the presence of outliers masks the structure of the rest of the data.
##' Outliers are defined by the 2.5 and 97.5 percentiles. 
##' @param method One of \code{PCA} (default), \code{MDS} or
##' \code{kpca}, defining if dimensionality reduction is done using
##' principal component analysis (see \code{\link{prcomp}}), classical
##' multidimensional scaling (see \code{\link{cmdscale}}) or kernel
##' PCA (see \code{kernlab::kpca}).
##' @param methargs A \code{list} of arguments to be passed when
##' \code{method} is called. If missing, the data will be scaled and
##' centred prior to PCA.
##' @param axsSwitch A \code{logical} indicating whether the axes should be
##' switched.
##' @param mirrorX A \code{logical} indicating whether the x axis should be mirrored? 
##' @param mirrorY A \code{logical} indicating whether the y axis should be mirrored? 
##' @param col A \code{character} of appropriate length defining colours.
##' @param pch A \code{character} of appropriate length defining point character.
##' @param cex Character expansion.
##' @param index A \code{logical} (default is FALSE), indicating if the feature
##' indices should be plotted on top of the symbols.
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
##' @seealso \code{\link{addLegend}} to add a legend to
##' \code{plot2D_v1} figures and \code{\link{plotDist}} for
##' alternative graphical representation of quantitative organelle
##' proteomics data.
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' pRoloc:::plot2D_v1(dunkley2006, fcol = NULL)
##' pRoloc:::plot2D_v1(dunkley2006, fcol = NULL, method = "kpca")
##' pRoloc:::plot2D_v1(dunkley2006, fcol = NULL, method = "kpca",
##'                    methargs = list(kpar = list(sigma = 1)))
##' pRoloc:::plot2D_v1(dunkley2006, fcol = "markers")
##' pRoloc:::addLegend_v1(dunkley2006,
##'                       fcol = "markers",
##'                       where = "topright",
##'                       cex = 0.5, bty = "n", ncol = 3)
##' title(main = "plot2D example")
plot2D_v1 <- function(object,
                      fcol = "markers",
                      fpch,
                      unknown = "unknown",
                      dims = 1:2,
                      alpha,
                      score = 1, ## TODO
                      outliers = TRUE,
                      method = c("PCA", "MDS", "kpca"),
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
    if (all(substr(stockcol,1 ,1) == "#") & !missing(alpha)) {
        if (alpha < 0)
            alpha <- 0
        if (alpha > 1)
            alpha <- 1
        alpha <- as.hexmode(as.integer(alpha * 255))
        stockcol <- paste0(stockcol, alpha)
    }  
    if (!is.null(fcol) && !fcol %in% fvarLabels(object))
        stop("'", fcol, "' not found in feature variables.")
    if (!missing(fpch) && !fpch %in% fvarLabels(object))
        stop("'", fpch, "' not found in feature variables.")
    method <- match.arg(method)
    if (length(dims) > 2) {
        warning("Using first two dimensions of ", dims)
        dims <- dims[1:2]
    }
    k <- max(dims)
    if (any(is.na(exprs(object)))) {
        narows <- unique(which(is.na(exprs(object)),
                               arr.ind = TRUE)[, "row"])
        object <- object[-narows, ]
        if (nrow(object) == 0)
            stop("No rows left after removing NAs!")
        else
            warning("Removed ", length(narows), " row(s) with 'NA' values.")    
    } 
    if (method == "PCA") {
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
    } else { ## kpca
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
    } 
    if (plot) {
        if (axsSwitch) {
            .data <- .data[, 2:1]
            .tmp <- .xlab
            .xlab <- .ylab
            .ylab <- .tmp
        }
        if (mirrorX)
            .data[, 1] <- -.data[, 1]
        if (mirrorY)
            .data[, 2] <- -.data[, 2]
        
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
            ukn <- fData(object)[, fcol] == unknown
        } else {
            ukn <- rep(TRUE, nrow(.data))
        }
        
        if (!is.null(fcol)) {
            .fcol <- factor(fData(object)[, fcol])    
            col <- stockcol[as.numeric(.fcol)]
            col[ukn] <- unknowncol
        }
        if (!missing(fpch)) {   
            .fpch <- factor(fData(object)[, fpch])
            pch <- stockpch[as.numeric(.fpch)]
        } else {
            pch <- rep(19, nrow(.data))
        }
        pch[ukn] <- unknownpch
        if (!outliers) {
            qntls <- apply(.data, 2, quantile, c(0.025, 0.975))
            selqtls <- .data[, 1] > qntls[1, 1] &
                .data[, 1] < qntls[2, 1] &
                    .data[, 2] > qntls[1, 2] &
                        .data[, 2] < qntls[2, 2] 
            .data <- .data[selqtls, ]
            ukn <- ukn[selqtls]
        }
        
        if (is.null(fcol)) {
            plot(.data, xlab = .xlab, ylab = .ylab)
        } else {
            plot(.data, xlab = .xlab, ylab = .ylab,
                 type = "n", ...)
            points(.data[ukn, ], col = col[ukn],
                   pch = pch[ukn], cex = cex[ukn], ...)
            points(.data[!ukn, ], col = col[!ukn],
                   pch = pch[!ukn], cex = cex[!ukn], ...)
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

##' Adds a legend to a \code{\link{plot2D_v1}} figure.
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
##' @param ... Additional parameters passed to \code{\link{legend}}.
##' @return Invisibly returns \code{NULL}
##' @author Laurent Gatto
addLegend_v1 <- function(object,
                         fcol = "markers",
                         where = "other",
                         col,
                         ...) {
  if (is.null(fcol))
    fcol <- "markers"
  if (!fcol %in% fvarLabels(object))
    stop("'", fcol, "' not found in feature variables.")
  if (missing(col)) {
    stockcol <- getStockcol()
  } else {
    stockcol <- col
  }  
  unknowncol <- getUnknowncol()
  unknownpch <- getUnknownpch()
  txt <- levels(factor(fData(object)[, fcol]))
  col <- stockcol[seq_len(length(txt))]
  col[txt == "unknown"] <- unknowncol
  pch <- rep(19, length(txt))
  pch[txt == "unknown"] <- unknownpch
  if (where == "other") {
    dev.new()
    plot(0, type = "n", bty = "n",
         xaxt = "n", yaxt = "n",
         xlab = "", ylab = "")
    legend("center",
           txt, col = col,
           pch = pch,
           bty = "n")
  } else if (where %in% c("bottomright",
                          "bottom", "bottomleft",
                          "left", "topleft", "top",
                          "topright", "right",
                          "center")) {
    legend(where, txt, col = col,
           pch = pch, ...)
  }
  invisible(NULL)
}


##' Produces a line plot showing the feature abundances
##' across the fractions.
##'
##' @title Plots the distribution of features across fractions
##' @param object An instance of class \code{MSnSet}.
##' @param markers A \code{character}, \code{numeric} or \code{logical}
##' of appropriate length and or content used to subset \code{object}
##' and define the organelle markers.
##' @param mcol A \code{character} define the colour of the marker features.
##' Default is \code{"steelblue"}. 
##' @param pcol A \code{character} define the colour of the non-markers features.
##' Default is \code{"grey90"}. 
##' @param alpha A numeric defining the alpha channel (transparency)
##' of the points, where \code{0 <= alpha <= 1}, 0 and 1 being completely
##' transparent and opaque.
##' @param lty Vector of line types for the marker profiles. Default
##' is 1 (solid). See \code{\link{par}} for details.
##' @param fractions An optional \code{character} defining the \code{phenoData}
##' variable to be used to label the fraction along the x axis. If missing, the
##' \code{phenoData} variables are searched for a match to \code{fraction}.
##' If no match is found, the fractions are labelled as numericals.
##' @param ylim A numeric vector of length 2, giving the y coordinates range.
##' @param ... Additional parameters passed to \code{\link{plot}}.
##' @return Used for its side effect of producing a feature distribution
##' plot. Invisibly returns \code{NULL}.
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
               pch = 1,
               col = mcol,
               type = "b",
               lty = lty,
               ...)
    }
  invisible(NULL)
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


##' Plot organelle assignment data and results.
##' 
##' Generate 2 dimensional or feature distribution plots to illustrate
##' localistation clusters. In \code{plot2D}, rows containing \code{NA}
##' values are removed prior to dimention reduction.
##'
##' Note that \code{plot2D} has been update in version 1.3.6 to
##' support more more orgnalle classes than colours defined in
##' \code{\link{getStockcol}}. In such cases, the defauly colours are
##' recycled using the default plotting characters defined in
##' \code{\link{getStockpch}}. See the example for an
##' illustration. The \code{alpha} argument is also depreciated in
##' version 1.3.6. Use \code{setStockcol} to set colours with
##' transparency instead. See example below.
##' 
##' @param object An instance of class \code{MSnSet}.
##' @param fcol Feature meta-data label (fData column name) defining
##' the groups to be differentiated using different colours. Default
##' is \code{markers}. Use \code{NULL} to suppress any colouring.
##' @param fpch Featre meta-data label (fData column name) desining
##' the groups to be differentiated using different point symbols.
##' @param unknown A \code{character} (default is \code{"unknown"})
##' defining how proteins of unknown localisation are labelled.
##' @param dims A \code{numeric} of length 2 defining the dimensions
##' to be plotted. Always 1:2 for MDS.
##' @param score A numeric specifying the minimum organelle assignment score
##' to consider features to be assigned an organelle. (not yet implemented).
##' @param method One of \code{"PCA"} (default), \code{"MDS"},
##' \code{"kpca"} or \code{"t-SNE"}, defining if dimensionality
##' reduction is done using principal component analysis (see
##' \code{\link{prcomp}}), classical multidimensional scaling (see
##' \code{\link{cmdscale}}), kernel ##' PCA (see \code{kernlab::kpca})
##' or t-SNE (see \code{tsne::tsne}). \code{"scree"} can also be used
##' to produce a scree plot.
##' @param methargs A \code{list} of arguments to be passed when
##' \code{method} is called. If missing, the data will be scaled and
##' centred prior to PCA.
##' @param axsSwitch A \code{logical} indicating whether the axes should be
##' switched.
##' @param mirrorX A \code{logical} indicating whether the x axis should be mirrored? 
##' @param mirrorY A \code{logical} indicating whether the y axis should be mirrored? 
##' @param col A \code{character} of appropriate length defining colours.
##' @param pch A \code{character} of appropriate length defining point character.
##' @param cex Character expansion.
##' @param index A \code{logical} (default is FALSE), indicating of the feature
##' indices should be plotted on top of the symbols.
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
##' figures and \code{\link{plotDist}} for alternative graphical
##' representation of quantitative organelle proteomics data.
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' plot2D(dunkley2006, fcol = NULL)
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
plot2D <- function(object,
                   fcol = "markers",
                   fpch,
                   unknown = "unknown",
                   dims = 1:2,
                   score = 1, ## TODO
                   method = c("PCA", "MDS", "kpca", "t-SNE", "scree"),
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

    if (!is.null(fcol) && !fcol %in% fvarLabels(object))
        stop("'", fcol, "' not found in feature variables.")
    if (!missing(fpch) && !fpch %in% fvarLabels(object))
        stop("'", fpch, "' not found in feature variables.")
    method <- match.arg(method)
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
        library("tsne")
        if (missing(methargs))
            .data <- tsne(exprs(object), k = k)
        else .data <- do.call(tsne,
                              c(list(X = exprs(object)),
                                k = k, 
                                methargs))
        .data <- .data[, dims]
        .xlab <- paste("Dimension 1")
        .ylab <- paste("Dimension 2")
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
    } else { ## kpca
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
    } 
    if (plot) {
        if (axsSwitch) {
            .data <- .data[, 2:1]
            .tmp <- .xlab
            .xlab <- .ylab
            .ylab <- .tmp
        }
        if (mirrorX)
            .data[, 1] <- -.data[, 1]
        if (mirrorY)
            .data[, 2] <- -.data[, 2]
        
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
##' @param ... Additional parameters passed to \code{\link{legend}}.
##' @return Invisibly returns \code{NULL}
##' @author Laurent Gatto
addLegend <- function(object,
                      fcol = "markers",
                      where = c("other", "bottomright", "bottom",
                          "bottomleft", "left", "topleft",
                          "top", "topright", "right", "center"),
                      col, 
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
    legend(where, txt, col = col, pch = pch, ...)    
    invisible(NULL)
}


##' Highlights a set of features of interest given as a
##' \code{FeaturesOfInterest} instance on a PCA plot produced by
##' code{plot2D}. If none of the features of interest are found in the
##' \code{MSnset}'s \code{featureNames}, an error is thrown.
##'
##' @title Highlight features of interest on a plot2D figure
##' @param object The main dataset described as an \code{MSnSet} or a
##' matrix with the coordinates of the features on the PCA plot
##' produced (and invisibly returned) by \code{plot2D}.
##' @param foi An instance of \code{\linkS4class{FeaturesOfInterest}}.
##' @param args A named list of arguments to be passed to
##' \code{plot2D} if the PCA coordinates are to be calculated. Ignored
##' if the PCA coordinates are passed directly, i.e. \code{object} is
##' a \code{matrix}.
##' @param ... Additional parameters passed to \code{points}.
##' @return NULL; used for its side effects.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data("tan2009r1")
##' x <- FeaturesOfInterest(description = "A test set of features of interest",
##'                         fnames = featureNames(tan2009r1)[1:10],
##'                         object = tan2009r1)
##' .pca <- plot2D(tan2009r1)
##' highlightOnPlot(.pca, x, col = "red")
##' highlightOnPlot(tan2009r1, x, col = "red", cex = 1.5)
##'
##' .pca <- plot2D(tan2009r1, dims = c(1, 3))
##' highlightOnPlot(.pca, x, pch = "+")
##' highlightOnPlot(tan2009r1, x, args = list(dims = c(1, 3)))
##'
##' plot2D(tan2009r1, mirrorX = TRUE)
##' highlightOnPlot(.pca, x, pch = "+", args = list(mirrorX = TRUE))
highlightOnPlot <- function(object, foi, args = list(), ...) {
    if (!fnamesIn(foi, object))
        stop("None of the features of interest are present in the data.")
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
    if (!is.null(args$mirrorX) && args$mirrorX)
        .pca[, 1] <- -.pca[, 1]
    if (!is.null(args$mirrorY) && args$mirrorY)
        .pca[, 2] <- -.pca[, 2]    
    points(.pca[sel, 1], .pca[sel, 2], ...)
}
