## setClass("ClustRegRes",
##          representation(algorithm = "character",
##                         parameters = "list",
##                         results = "matrix",
##                         criteria = "character",
##                         algoparams = "list"))

## setMethod("show",
##           signature(object = "ClustRegRes"),
##           function(object) {
##             cat("Object of class \"", class(object), "\"\n", sep = "")
##             cat(" Algorithm:", object@algorithm, "\n")
##             cat(" Criteria:", paste(object@criteria, collapse = " "), "\n")
##             cat(" Parameters:\n")
##             for (i in 1:length(object@parameters)) {
##               cat(" ", names(object@parameters)[i],": ")
##               .param <- object@parameters[[i]]
##               .n <- length(.param)
##               if (.n < 7) cat(paste0(paste(.param, collapse = " "), "\n"))
##               else cat(paste0(paste(.param[1:2], collapse = " "),
##                               " ... ",
##                               paste(.param[(.n-1):.n], collapse = " "),
##                               "\n"))
##             }
##             invisible(NULL) 
##           })

## setMethod("getParams", "ClustRegRes",
##           function(object, criterion) {
##             if (missing(criterion)) {
##               criterion <- object@criteria[1]
##             } else {
##               criterion <- match.arg(criterion,
##                                      choices = object@criteria,
##                                      several.ok = FALSE)
##             }
##             idx <- which.min(object@results[, criterion])
##             return(object@results[idx, names(object@parameters)])
##           })

## setMethod("plot", c("ClustRegRes", "missing"),
##           function(x, y, ...) {
##             cls <- c("steelblue", "red")
##             res <- x@results[, x@criteria]            
##             matplot(res, pch = 1, type = "b",
##                     xlab = "k", ylab = "IC",
##                     col = cls)
##             grid()
##             xx <- apply(res, 2, which.min)
##             yy <- apply(res, 2, min)
##             points(xx, yy, col = cls, pch = 19)            
##             legend("topright", c("BIC", "AIC"),
##                    pch = 19, col = cls, bty = "n")
##             invisible(NULL)            
##           })


## setMethod("levelPlot", "ClustRegRes", 
##           function(object, ...) {
##             res <- object@results[, object@criteria]
##             mins <- apply(res, 2, min)
##             for (i in 1:ncol(res))
##               res[, i] <- res[, i]/mins[i]
##             p <- levelplot(res, xlab = "k",
##                            ylab = NULL,
##                            ...)
##             p
##           })









