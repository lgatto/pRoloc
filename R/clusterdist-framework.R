######################################################################
## ClustDist: this class summarises the algorithm information from
## running the clustDist algorithm. Information such as (i) the number 
## of k's tested for the kmeans, (ii) the mean and (iii) normalised, 
## pairwise Euclidean distances and (iv) cluster size per numer of 
## component clusters tested, for each GO term tested.
setClass("ClustDist",
         representation(k = "numeric",
                        dist = "list",
                        term = "character",
                        id = "character",
                        nrow = "numeric",
                        clustsz = "list",
                        components = "vector",
                        fcol = "character"))

setMethod("show",
          signature(object = "ClustDist"),
          function(object) {
            cat("Object of class \"",class(object),"\"\n",sep="")
            cat("fcol = ", object@fcol, "\n")
            cat(" term = ", object@term, "\n")
            cat(" id = ", object@id, "\n")
            cat(" nrow = ", object@nrow, "\n")
            cat("k's tested:", object@k, "\n")
            for (i in 1:length(object@k))
              cat("  Size: ", 
                  paste(object@clustsz[[i]], collapse = ", ")
                  , "\n")
            dfr <- getClusterInfo(object)
            cat("Clusters info:\n")
            print(dfr)
            invisible(NULL)
          })

##      x =  Object of "ClustDist"
##      y =  Object of "MSnSet"
## method =  One of "norm" or "mean", the default is "norm", indicating whether
##           the mean distance, or mean normalised distance per k clusters 
##           should be calculated
##      p =  Normalisation factor. Default is 1/3.
##      nchar = Maximum number of characters of GO ID, before their truncation. Default is 40.
setMethod("plot", c("ClustDist", "MSnSet"),
          function(x, y,
                   method = "norm",
                   p = 1/3,
                   nchar = 40
          ) {
            opar <- par(no.readonly = TRUE)
            on.exit(par(opar))
            if (!(method == "norm" | method == "mean"))
              stop("method must be one of 'norm' or 'mean' see ?plotClustDist for details")
            if (method == "norm")
              clusterdists <- normDist(x, best = FALSE, p = p)
            if (method == "mean")
              clusterdists <- meanDist(x, best = FALSE) ## CHECK HERE
            clusterdists <- signif(clusterdists, 3)
            numk <- length(x@k)
            if (numk > 3)
              par(mfrow = c(floor(sqrt(numk)), ceiling(sqrt(numk))),
                  oma = c(0, 0, 2, 0))
            else
              par(mfrow = c(1, numk), oma = c(0, 0, 2, 0))
            title <- x@id
            for (i in 1:length(x@k)) {
              fData(y)$.tmp <- "unknown"
              components <- x@components[[i]]
              fData(y)[names(components), ".tmp"] <- as.character(components)
              plot2D(y, ".tmp", main = paste("distance = ", clusterdists[i]))
            }
            mtext(title, outer = TRUE, cex = 1, font = 2)
          })

###################################################################
## ClustDistList: container for multiple ClustDist objects.
.ClustDistList <-
  setClass("ClustDistList",
           slots = c(x = "list",
                     log = "list"),
           contains = "Versioned",
           prototype = prototype(
             new("Versioned",
                 versions = c(ClustDistList = "0.1.0"))),
           validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (!listOf(object@x, "ClustDist"))
               msg <- validMsg(msg, "Not all items are ClustDist objects")
             nvals <- sapply(object@x, validObject)
             if (!all(nvals))
               msg <- validMsg(msg,
                               paste(sum(!nvals),
                                     "ClustDists are not valid."))
             if (is.null(msg)) TRUE
             else msg
           })

ClustDistList <-
  function(x, log = list(call = match.call()))
    .ClustDistList(x = x, log = log)

setMethod("show", "ClustDistList",
          function(object) {
            cat("Instance of class '", class(object), "' containig ",
                length(object), " objects.\n", sep = "")
          })

setMethod("length", "ClustDistList", function(x) length(x@x))

setMethod("names", "ClustDistList", function(x) names(x@x))

setMethod("[", c("ClustDistList", "ANY", "missing", "missing"),
          function(x, i, j = "missing", drop = "missing")
            .ClustDistList(x = msnsets(x)[i]))

setMethod("[[", c("ClustDistList", "ANY", "missing"),
          function(x, i, j = "missing", drop = "missing") {
            if (length(i) != 1)
              stop("subscript out of bounds")
            msnsets(x)[[i]]
          })

setReplaceMethod("names", "ClustDistList",
                 function(x, value) {
                   names(x@x) <- value
                   x
                 })

clustdists <- function(object) object@x

setMethod("lapply", "ClustDistList",
          function(X, FUN, ...) {
            ans <- lapply(clustdists(X), FUN, ...)
            if (listOf(ans, "ClustDist"))
              ans <- ClustDistList(ans)
            ans
          })

setMethod("sapply", "ClustDistList",
          function(X, FUN, ...) {
            ans <- sapply(clustdists(X), FUN, ...)
            if (listOf(ans, "ClustDist"))
              ans <- ClustDistList(ans)
            ans
          })

##      x =  Object of "ClustDistList"
## method =  One of "norm" or "mean", the default is "norm", indicating whether
##           the mean distance, or mean normalised distance per k clusters 
##           should be calculated
##      p =  Normalisation factor. Default is 1/3.
##  nchar = Maximum number of characters of GO ID, before their truncation. Default is 40.
##   ...  = Arguments passed to axis
setMethod("plot", c("ClustDistList", "missing"),
          function(x, y,
                   method = "norm",
                   p = 1/3,
                   main = "",
                   mai = c(5, 1.1, .5, .56),
                   mar = c(15, 6, 2, 2),
                   nchar = 40,
                   ...
          ) {
            opar <- par(no.readonly = TRUE)
            on.exit(par(opar))
            if (!(method == "norm" | method == "mean"))
              stop("method must be one of 'norm' or 'mean' see ?ClustDistList for details")
            if (method == "norm")
              clusterdists <- lapply(x, normDist, best = FALSE, p = p)
            if (method == "mean")
              clusterdists <- lapply(x, meanDist, best = FALSE)
            orgnames <- sapply(x, function(z) z@id)
            ll <- sapply(orgnames, nchar)
            to.shorten <- which(ll > nchar)
            if (length(to.shorten) > 0) {
              replace.with <- paste0(substr(orgnames[to.shorten], 1, 37), "...")
              orgnames[to.shorten] <- replace.with
            }
            min.dist <- sapply(clusterdists, min, na.rm = TRUE)
            oo <- order(min.dist, decreasing = TRUE)
            names(clusterdists) <- orgnames
            orgnames <- orgnames[oo]
            clusterdists <- clusterdists[oo]
            if (method == "norm")
              ylab <- "Mean normalised \nEuclidean distance (per k)"
            else
              ylab <- "Mean Euclidean \ndistance (per k)"
            par(mai = c(5, 1.1, .5, .56), mar = c(15, 6, 2, 2) + 0.1)
            boxplot(clusterdists, type = "p", pch = 19, cex = .5,
                    xaxt = "n", yaxt = "n", axt = "n", xlab = NA,
                    ylab = ylab, main = main, cex.main = .8)
            axis(1, at = 1:length(oo), labels = orgnames,
                 las = 2, cex.axis = .7, ...)
            axis(2, at = seq(0, max(unlist(clusterdists), na.rm = TRUE), by = .05),
                 las = 2, cex.axis = .9)
          })
