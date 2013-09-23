##' Plot organelle assignment data and results.
##' 
##' Generate 2 dimensional or feature distribution plots to illustrate
##' localistation clusters. In \code{plot2D}, rows containing \code{NA}
##' values are removed prior to dimention reduction.
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
##' to be plotted, i.e the PC/MDS axes. 
##' @param alpha A numeric defining the alpha channel (transparency)
##' of the points, where \code{0 <= alpha <= 1}, 0 and 1 being completely
##' transparent and opaque.
##' @param score A numeric specifying the minimum organelle assignment score
##' to consider features to be assigned an organelle. (not yet implemented).
##' @param outliers A logical specifying whether outliers should be plotted
##' or ignored (default is TRUE, i.e. all points are plotted). Useful when 
##' the presence of outliers masks the structure of the rest of the data.
##' Outliers are defined by the 2.5 and 97.5 percentiles. 
##' @param method One of \code{PCA} (default) or \code{MDS}, defining
##' if dimensionality reduction is done using principal component
##' analysis (see \code{\link{prcomp}}) or classical multidimensional
##' scaling  (see \code{\link{cmdscale}}). 
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
##' plot2D(dunkley2006, fcol = NULL, method = "MDS")
##' plot2D(dunkley2006, fcol = "markers")
##' addLegend(dunkley2006,
##'           fcol = "markers",
##'           where = "topright",
##'           cex = 0.5, bty = "n", ncol = 3)
##' title(main = "plot2D example")
plot2D <- function(object,
                   fcol = "markers",
                   fpch,
                   unknown = "unknown",
                   dims = 1:2,
                   alpha,
                   score = 1, ## TODO
                   outliers = TRUE,
                   method = c("PCA", "MDS"),
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
    .pca <- prcomp(exprs(object), scale=TRUE, center=TRUE)
    .data <- .pca$x[, dims]
    .vars <- (.pca$sdev)^2
    .vars <- (.vars / sum(.vars))[dims]
    .vars <- round(100 * .vars, 2)
    .xlab <- paste0("PC", dims[1], " (", .vars[1], "%)")
    .ylab <- paste0("PC", dims[2], " (", .vars[2], "%)")
  } else { ## MDS
    .data <- cmdscale(dist(exprs(object), 
                           method = "euclidean",
                           diag = FALSE,
                           upper = FALSE,
                           p = 2),
                      eig = TRUE,
                      k = k)$points
    .data <- .data[, dims]
    .xlab <- paste("Dimension", dims[1])
    .ylab <- paste("Dimension", dims[2])    
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
        qntls <- apply(.data, 2, quantile, c(0.025, 0.971))
        selqtls <- .data[, 1] > qntls[1, 1] &
            .data[, 1] < qntls[2, 1] &
                .data[, 2] > qntls[1, 2] &
                    .data[, 2] < qntls[2, 2] 
        .data <- .data[selqtls, ]
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

##' Adds a legend to a \code{\link{plot2D}} figure.
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
##' @param fractions An optional \code{character} defining the \code{phenoData}
##' variable to be used to label the fraction along the x axis. If missing, the
##' \code{phenoData} variables are searched for a match to \code{fraction}.
##' If no match is found, the fractions are labelled as numericals.
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
##'          markers = featureNames(tan2009r1)[j],
##'          main = "Mitochondrion")
plotDist <- function(object,
                     markers,
                     mcol = "steelblue",                     
                     pcol = "grey90",
                     alpha = 0.3,
                     fractions,
                     ...) {
  .data <- exprs(object)
  .ylim <- range(.data)
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
  plot(0, ylim = .ylim,
       xlim = c(1, m),
       ylab = "Intensity",
       xlab = "Fractions",
       type = "n", xaxt = "n",
       ...)
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
               lty = "solid")
    }
  invisible(NULL)
}

