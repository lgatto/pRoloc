##' Plot organelle assignment data and results.
##' 
##' Generate 2 dimensional or feature distribution plots to illustrate
##' localistation clusters. In \code{plot2D}, rows containing \code{NA}
##' values are removed prior to dimention reduction.
##' 
##' TODO
##' 
##' @aliases plot2D plotDist addLegend
##' @param object An instance of class \code{MSnSet}.
##' @param fcol
##' @param fpch
##' @param unknown
##' @param dims
##' @param alpha
##' @param score
##' @param method
##' @param axsSwitch
##' @param mirrorX
##' @param mirrorY
##' @param col
##' @param pch
##' @param where (addLegend only)
##' @param markers
##' @param mcol
##' @param pcol
##' @param fractions
##' @param \dots
##' @return Used for their side effects of generating plots and
##' adding legends. Invisibly returns the 2d data.
##' @author Laurent Gatto <lg390@@cam.ac.uk>
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' plot2D(dunkley2006, fcol = NULL)
##' plot2D(dunkley2006, fcol = "markers")
##' addLegend(dunkley2006,
##'           fcol = "markers",
##'           where = "topright",
##'           cex = 0.5, bty = "n", ncol = 3)
##' title(main = "plot2D example")
##' data(tan2009r1)
##' j <- which(fData(tan2009r1)$markers == "mitochondrion")
##' i <- which(fData(tan2009r1)$PLSDA == "mitochondrion")
##' plotDist(tan2009r1[i, ],
##'          markers = featureNames(tan2009r1)[j],
##'          main = "Mitochondrion")
##' 
plot2D <- function(object,
                   fcol = "markers",
                   fpch,
                   unknown = "unknown",
                   dims = 1:2,
                   alpha,
                   score = 1, ## TODO
                   method = c("PCA", "MDS"),
                   axsSwitch = FALSE,
                   mirrorX = FALSE,
                   mirrorY = FALSE,
                   col,
                   pch,
                   cex,
                   identify = FALSE,
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
    .data <- cmdscale(dist(as.matrix(exprs(object)), 
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
  if (identify) {
    ids <- identify(.data[, 1], .data[, 2],
                    rownames(.data))
    return(ids)    
  }
  invisible(.data)
}


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
}


plotDist <- function(object,
                     markers,
                     mcol = "steelblue",                     
                     proteins,
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
  axis(1, at = seq_len(m), label = .frac)
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
}

