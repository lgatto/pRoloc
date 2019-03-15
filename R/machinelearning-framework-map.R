## Documented with tagmMapTrain and tagmMapPredict in
## machinelearning-functions-tagm-map.R

##' @slot method A `character()` storing the TAGM method name.
##' @slot priors A `list()` with the priors for the parameters
##' @slot seed An `integer()` with the random number generation seed.
##' @slot posteriors A `list()` with the updated posterior parameters
##'     and log-posterior of the model.
##' @slot datasize A `list()` with details about size of data
##' @md
##' @aliases class:MAPParams MAPParams-class MAPParams
##' @rdname tagm-map
##' @author Laurent Gatto
setClass("MAPParams",
         representation(method = "character",
                        priors = "list",
                        seed = "integer",
                        posteriors = "list",
                        datasize = "list"))

##' @rdname tagm-map
setMethod("show", "MAPParams",
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep = "")
            cat(" Method:", object@method, "\n")
            invisible(NULL)
          })

##' The `logPosteriors` function can be used to extract the log-posteriors at
##' each iteration of the EM algorithm to check for convergence.
##' @param x An object of class `MAPParams`.
##' @rdname tagm-map
logPosteriors <- function(x) {
    stopifnot(inherits(x, "MAPParams"))
    x@posteriors$logposterior
}