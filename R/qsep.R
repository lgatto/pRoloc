.QSep <- setClass("QSep",
                  contains = c("Versioned"),
                  slots = c(x = "matrix",
                            xnorm = "matrix",
                            object = "character"),
                  prototype = prototype(
                      new("Versioned",
                          versions = c(QSep = "0.1.0"))),
                  validity = function(object) {
                      msg <- validMsg(NULL, NULL)
                      if (!identical(dim(object@x), dim(object@xnorm)))
                          msg <- validMsg(msg, "Dimensions don't match.")
                      if (!identical(dimnames(object@x), dimnames(object@xnorm)))
                          msg <- validMsg(msg, "Dimension names don't match.")
                      if (is.null(msg)) TRUE
                      else msg
                  })

QSep <- function(object) {
    objname <- MSnbase:::getVariableName(match.call(), "object")
    ## only consider markers
    hlm <- markerMSnSet(object)
    ## vector of markers
    um <- unique(getMarkers(hlm, verbose = FALSE))
    ## answer is a square matrix
    ans <- diag(length(um))
    colnames(ans) <- rownames(ans) <- um
    ## euclidean distance between all markers
    hlmd <- dist(exprs(hlm))
    hlmd <- as.matrix(hlmd)
    diag(hlmd) <- NA
    ## mean distance between all pairs of subcellular marker clusters
    tmp <- apply(expand.grid(um, um), 1,
                 function(.um) {
                     sel1 <- fData(hlm)$markers == .um[1]
                     sel2 <- fData(hlm)$markers == .um[2]
                     ans[.um[1], .um[2]] <<- mean(hlmd[sel1, sel2],
                                                  na.rm = TRUE)
                 })
    res <- .QSep(x = ans,
                 xnorm = ans / diag(ans),
                 object = objname)
    if (validObject(res)) res
}

qsep <- function(object, norm = TRUE) {
    stopifnot(inherits(object, "QSep"))
    if (norm) object@xnorm
    else object@x
}

setMethod("names", "QSep",
          function(x) colnames(x@x))

setReplaceMethod("names", signature(x = "QSep", value = "character"),
                 function(x, value) {
                     if (length(value) != nrow(x@x))
                         stop("New names must be of length ", nrow(x@x))
                     colnames(x@x) <- colnames(x@xnorm) <-
                         rownames(x@x) <- rownames(x@xnorm) <-
                         value
                     x
                 })

setMethod("show", "QSep",
          function(object) {
              cat("Object of class '", class(object), "'.\n", sep = "")
              cat(" Data:", object@object, "\n")
              cat(" With", nrow(object@x), "sub-cellular clusters.\n")
              invisible(NULL)
          })

setMethod("summary", "QSep",
          function(object, ..., verbose = TRUE) {
              ans <- object@xnorm
              diag(ans) <- NA
              ans <- na.omit(as.numeric(ans))
              if (verbose)
                  print(summary(ans))
              invisible(ans)
          })

setMethod("levelPlot", "QSep",
          function(object) {
              pal <- colorRampPalette(c("blue", "white", "red"))
              myPanel <- function(x, y, z, ...) {
                  panel.levelplot(x, y, z, ...)
                  panel.text(x, y, ifelse(is.na(z), "", round(z, 1)))
              }
              levelplot(qsep(object),
                        col.regions = pal(50),
                        panel = myPanel,
                        xlab = "", ylab = "",
                        scales = list(x = list(cex = .8, rot = 45),
                                      y = list(cex = .8)))
          })

.plotQSep <- function(obj, ...) {
    args <- pairlist(...)
    x <- qsep(obj, norm = TRUE)
    n <- nrow(x)
    opar <- par(las = 1)
    on.exit(par(opar))
    if (!"ylim" %in% names(args))
        boxplot(x, horizontal = TRUE,
                col = "#00000010",
                ylim = c(min(x), max(x)), ...)
    else
        boxplot(x, horizontal = TRUE,
                col = "#00000010",
                ...)
    points(as.numeric(x),
           rep(1:n, each = n),
           col = "#00000090",
           pch = 19)
    points(diag(x), 1:n, col = "red", pch = 1, lwd = 2)
}

setMethod("plot", c("QSep", "missing"),
          function(x, y, ...)
              .plotQSep(x, ...))
