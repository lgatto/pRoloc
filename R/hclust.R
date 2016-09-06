##' This functions calculates an average protein profile for each
##' marker class (proteins of unknown localisation are ignored) and
##' then generates a dendrogram representing the relation between
##' marker classes. The colours used for the dendrogram labels are
##' taken from the default colours (see \code{\link{getStockcol}}) so
##' as to match the colours with other spatial proteomics
##' visualisations such as \code{\link{plot2D}}.
##'
##' @title Draw a dendrogram of subcellular clusters
##' @param object An instance of class \code{MSnSet}.
##' @param fcol Feature meta-data label (fData column name) defining
##'     the groups to be differentiated using different
##'     colours. Default is \code{markers}.
##' @param distargs A \code{list} of arguments to be passed to the
##'     \code{\link[stats]{dist}} function.
##' @param hclustargs A \code{list} of arguments to be passed to the
##'     \code{\link[stats]{hclust}} function.
##' @param plot A \code{logical} defining whether the dendrogram
##'     should be plotted. Default is \code{TRUE}.
##' @param ... Additional parameters passed to
##'     \code{\link[stats]{stats::plot.dendrogram}}.
##' @return Invisibly returns a matrix of average occupancy profiles
##'     for all marker classes defined in \code{fcol}. 
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' mrkHClust(dunkley2006)
mrkHClust <- function(object, fcol = "markers",
                      distargs, hclustargs,
                      plot = TRUE,
                      ...) {
    object <- markerMSnSet(object)
    mm <- getMarkers(object, fcol = fcol, verbose = FALSE)
    umm <- levels(factor(mm))
    fmat <- matrix(NA, nrow = length(umm), ncol = ncol(object))
    rownames(fmat) <- umm
    colnames(fmat) <- sampleNames(object)
    for (m in umm)
        fmat[m, ] <- colMeans(exprs(object[mm == m, ]), na.rm = TRUE)
    ## handling additional dist args
    if (missing(distargs)) distargs <- list(x = fmat)
    else distargs <- c(list(x = fmat), distargs)
    d <- do.call(dist, distargs)
    ## handling additional hclust args
    if (missing(hclustargs)) hclustargs <- list(d = d)
    else hclustargs <- c(list(d = d), hclustargs)
    hc <- do.call(hclust, hclustargs)
    hc <- as.dendrogram(hc)
    i <- match(labels(hc), umm)
    hc <- set(hc, "labels_colors", getStockcol()[i])
    if (plot) plot(hc, ...)
    invisible(fmat)
}
