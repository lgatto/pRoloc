## This will most likely be moved to MSnbase

## library("pRoloc")
## library("pRolocdata")
## data(tan2009r1)
## library("digest")

.FeaturesOfInterest <-
    setClass("FeaturesOfInterest",
         slots = c(description = "character",
             objpar = "list",
             fnames = "character",
             date = "character"))

setGeneric("foi",
           function(object, ...) standardGeneric("foi"))

setMethod("foi", "FeaturesOfInterest",
          function(object, ...) object@fnames)

setMethod("description", "FeaturesOfInterest",
          function(object, ...) object@description)

setMethod("length", "FeaturesOfInterest",
         function(x) length(x@fnames))

setMethod("show", "FeaturesOfInterest",
          function(object) {
              cat("Object of class \"", class(object), "\"\n", sep="")
              cat(" Created on", object@date, "\n")
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


FeaturesOfInterest <- function(fnames,
                            description,
                            object) {
    fns <- featureNames(object)
    if (!all(fnames %in% fns)) {
        fx <- sum(!fnames %in% fns)
        stop(fx, " feature(s) of interest absent from your object's feature names:\n   ",
             paste(fnames[fx], collapse = ", "))
    }
    .FeaturesOfInterest(description = description,
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


setGeneric("addFeaturesOfInterest",
           function(x, y) standardGeneric("addFeaturesOfInterest"))

setMethod("addFeaturesOfInterest",
          c("FeaturesOfInterest", "FoICollection"),
          function(x, y) {              
              y@foic <- c(y@foic, x)
              return(y)
          })


setGeneric("rmFeaturesOfInterest",
           function(object, i) standardGeneric("rmFeaturesOfInterest"))

setMethod("rmFeaturesOfInterest",
          c("FoICollection", "numeric"),
          function(object, i) {
              object@foic <- object@foic[-i]
              return(object)
          })


setMethod("description", "FoICollection",
          function(object, ...) sapply(foi(object), description))


setGeneric("fromIdentical",
           function(x, y, ...) standardGeneric("fromIdentical"))

setMethod("fromIdentical",
          c("FeaturesOfInterest", "FeaturesOfInterest"),
          function(x, y) x@objpar$digest == y@objpar$digest)

setMethod("fromIdentical",
          c("FoICollection", "missing"),
          function(x) {
              dgsts <- sapply(foi(x), function(xx) xx@objpar$digest)
              length(unique(dgsts)) == 1
          })
              

setGeneric("fromEqual",
           function(x, y, ...) standardGeneric("fromEqual"))


setMethod("fromEqual",
          c("FeaturesOfInterest", "FeaturesOfInterest"),
          function(x, y)
          x@objpar$nrow == y@objpar$nrow & x@objpar$ncol == y@objpar$ncol)

setMethod("fromEqual",
          c("FoICollection", "missing"),
          function(x) {
              nc <- sapply(foi(x), function(xx) xx@objpar$ncol)
              nr <- sapply(foi(x), function(xx) xx@objpar$nrow)
              length(unique(nc)) == 1 & length(unique(nr)) == 1
          })



## x <- FeaturesOfInterest(description = "A test set of features of interest",
##                      fnames = featureNames(tan2009r1)[1:10],
##                      object = tan2009r1)

## y <- FeaturesOfInterest(description = "A second test set of features of interest",
##                        fnames = featureNames(tan2009r1)[111:113],
##                        object = tan2009r1)

## xx <- FoICollection()

## xx <- addFeaturesOfInterest(x, xx)
## xx <- addFeaturesOfInterest(y, xx)


## NOTES:
##     on application of a foi on an MSnSet:
##     if .@objpar$name differs -> message
##     if .@objpar$ncol differs -> warning
##     if .@objpar$nrow differs -> warning
##     if any of !foi(.) %in% featureNames(.) -> warning
## if no overlap -> error


highlightOnPlot <- function(object, foi, ...) {
    ## TODO: checks above
    .pca <- plot2D(object, plot = FALSE)
    sel <- featureNames(object) %in% foi(foi)
    points(.pca[sel, 1], .pca[sel, 2], ...)
}

## plot2D(tan2009r1)
## highlightOnPlot(tan2009r1, x, col = "red")
## highlightOnPlot(tan2009r1, y, col = "green", pch = "+")

