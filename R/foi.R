.FeatsOfInterest <-
    setClass("FeatsOfInterest",
         slots = c(description = "character",
             objpar = "list",
             fnames = "character",
             date = "character"))

setGeneric("foi",
           function(object, ...) standardGeneric("foi"))

setMethod("foi", "FeatsOfInterest",
          function(object, ...) object@fnames)

setMethod("description", "FeatsOfInterest",
          function(object, ...) object@description)

setMethod("length", "FeatsOfInterest",
         function(x) length(x@fnames))

setMethod("show", "FeatsOfInterest",
          function(object) {
              cat("Object of class \"", class(object), "\"\n", sep="")
              cat(" Description:\n")
              cat(strwrap(object@description,
                          width = 0.65 * getOption("width"),
                          indent = 2, exdent = 2), sep = "\n")
              n <- length(object)
              cat(" ", n, " features of interest:\n", sep = "")
              if (n > 5) {
                  cat("  ", paste(object@fnames[1:2], collapse = ", "),
                      " ... ", paste(object@fnames[(n-1):n], collapse = ", "))
                  } else {
                      cat("  ", paste(object@fnames, collapse = ", "))
                  }
              cat("\n")
          })


FeatsOfInterest <- function(fnames,
                            description,
                            object) {
    fns <- featureNames(object)
    if (!all(fnames %in% fns)) {
        fx <- sum(!fnames %in% fns)
        stop(fx, " feature(s) of interest absent from your object's feature names:\n   ",
             paste(fnames[fx], collapse = ", "))
    }
    .FeatsOfInterest(description = description,
                     fnames = fnames,
                     date = date(),
                     objpar = list(
                         ncol = ncol(object),
                         nrow = nrow(object),
                         name = MSnbase:::getVariableName(
                             match.call(),
                             "object"),
                         digest = digest(object)))
}

.FoICollection <- setClass("FoICollection",
                           slots = c(foic = "list"))

setGeneric("FoICollection",
           function(object, ...) standardGeneric("FoICollection"))

setMethod("FoICollection", "missing",
          function(object, ...)
          .FoICollection(foic = list()))

setMethod("FoICollection", "list",
          function(object, ...) .FoICollection(foic = object))

setMethod("length", "FoICollection",
         function(x) length(x@foic))

setMethod("show", "FoICollection",
          function(object)
          cat("A collection of ", length(object),
              " features of interest.\n", sep = ""))

setMethod("foi", "FoICollection",
          function(object, n) {
              if (missing(n)) return(object@foic)
              else return(object@foic[n])
          })


setGeneric("addFeatsOfInterest",
           function(x, y) standardGeneric("addFeatsOfInterest"))

setMethod("addFeatsOfInterest",
          c("FeatsOfInterest", "FoICollection"),
          function(x, y) {              
              y@foic <- c(y@foic, x)
              return(y)
          })


setGeneric("rmFeatsOfInterest",
           function(object, i) standardGeneric("rmFeatsOfInterest"))

setMethod("rmFeatsOfInterest",
          c("FoICollection", "numeric"),
          function(object, i) {
              object@foic <- object@foic[-i]
              return(object)
          })


setMethod("description", "FoICollection",
          function(object, ...) sapply(foi(object), description))


## x <- FeatsOfInterest(description = "A test set of features of interest",
##                        fnames = featureNames(tan2009r1)[1:10],
##                        object = tan2009r1)
## y <- FeatsOfInterest(description = "A second test set of features of interest",
##                        fnames = featureNames(tan2009r1)[111:113],
##                        object = tan2009r1)

## xx <- FoICollection()
## xx <- addFeatsOfInterest(x, xx)
## xx <- addFeatsOfInterest(y, xx)

## NOTES:
##     on application of a foi on an MSnSet:
##     if .@objpar$name differs -> message
##     if .@objpar$ncol differs -> warning
##     if .@objpar$nrow differs -> warning
##     if any of !foi(.) %in% featureNames(.) -> warning

highlightOnPlot <- function(object, foi, ...) {
    ## TODO: checks above
    .pca <- plot2D(object, plot = FALSE)
    sel <- featureNames(object) %in% foi(foi)
    points(.pca[sel, 1], .pca[sel, 2], ...)
}

## plot2D(tan2009r1)
## highlightOnPlot(tan2009r1, x, col = "red")
## highlightOnPlot(tan2009r1, y, col = "green", pch = "+")

