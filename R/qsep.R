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

.meanMarkerDist <- function(x, markers) {
  um <- unique(markers)
  n <- length(um)
  sel <- lapply(um, "==", markers)
  ans <- matrix(NA_real_, nrow=n, ncol=n, dimnames=list(um, um))

  for (i in 1L:n) {
    for (j in 1L:n) {
      ans[i, j] <- mean(x[sel[[i]], sel[[j]]], na.rm=TRUE)
    }
  }
  ans
}

QSep <- function(object, fcol = "markers") {
    objname <- MSnbase:::getVariableName(match.call(), "object")
    ## only consider markers
    mobj <- markerMSnSet(object, fcol = fcol)
    ## vector of markers
    markers <- getMarkers(mobj, fcol = fcol, verbose = FALSE)
    ## euclidean distance between all markers
    mrkdist <- dist(exprs(mobj))
    mrkdist <- as.matrix(mrkdist)
    diag(mrkdist) <- NA

    ans <- .meanMarkerDist(mrkdist, markers)

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
          function(object, norm = TRUE, ...) {
              pal <- colorRampPalette(c("blue", "white", "red"))
              myPanel <- function(x, y, z, ...) {
                  panel.levelplot(x, y, z, ...)
                  panel.text(x, y, ifelse(is.na(z), "", round(z, 2)))
              }
              levelplot(t(qsep(object, norm = norm)),
                        col.regions = pal(50),
                        panel = myPanel,
                        xlab = "Reference cluster", ylab = "",
                        scales = list(x = list(cex = .8, rot = 45),
                                      y = list(cex = .8)),
                        ...)
          })

.plotQSep <- function(obj, norm = TRUE, ...) {
    args <- pairlist(...)
    x <- qsep(obj, norm = norm)
    n <- nrow(x)

    ## add enough space on the left for ytick labels
    rni <- max(strwidth(rownames(x), "inch"), na.rm = TRUE)
    mai <- par("mai")
    mai[2] <- mai[4] + rni + 0.1

    opar <- par(las = 1,
                mai = mai)
    on.exit(par(opar))
    if (!"ylim" %in% names(args))
        boxplot(x, horizontal = TRUE,
                col = "#00000010",
                ylim = c(min(x, na.rm = TRUE),
                         max(x, na.rm = TRUE)),
                ...)
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
          function(x, y, norm = TRUE, ...)
              .plotQSep(x, norm,
                        ...))
