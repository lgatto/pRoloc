##' @slot method `character(1)` describing the method.
##' @slot chains Object of class `MCMCChains` containing the full MCMC
##'     chain results.
##' @slot summary Object of class `MCMCSummary` the summarised MCMC
##'     results.
.MCMCParams <- setClass("MCMCParams",
                        slots = c(method = "character",
                                  chains = "MCMCChains",
                                  summary = "MCMCSummary"))

chains <- function(x) {
    stopifnot(inherits(x, "MCMCParams"))
    x@chains
}


setMethod("show", "MCMCParams",
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep = "")
            cat("Method:", object@method, "\n")
            cat("Number of chains:", length(object@chains), "\n")
            invisible(NULL)
          })


.MCMCSummary <- setClass("MCMCSummary",
                         slots = c(summary = "list"))

##' @slot K `integer(1)` indicating the number of components.
##' @slot D = "integer", ## number of samples.
##' @slot method `character(1)` defining the method used. Currently
##'     `"TAGM.MCMC"`, later also `"TAGM.GP"`.
##' @slot mk `matrix(K, D)`
##' @slot lambdak `numeric(K)`
##' @slot nuk `numeric(K)`
##' @slot sk `array(K, D, D)`
.ComponentParam <- setClass("ComponentParam",
                            slots = c(K = "integer",
                                      D = "integer",
                                      method = "character",
                                      mk = "matrix",
                                      lambdak = "numeric",
                                      nuk = "numeric",
                                      sk = "array"),
                            prototype = prototype(
                                method = "TAGM.MCMC"
                                ),
                            validity = function(object){
                                msg <- validMsg(NULL, NULL)
                                K <- object@K
                                D <- object@D
                                if (object@method != "TAGM.MCMC")
                                    msg <- validMsg(msg, "Wrong method")
                                if (!identical(dim(object@mk), c(K, D)))
                                    msg <- validMsg(msg, "Wrong dimensions: mk")
                                if (!identical(length(object@lambdak), K))
                                    msg <- validMsg(msg, "Wrong length: lambdak")
                                if (!identical(length(object@nuk), K))
                                    msg <- validMsg(msg, "Wrong length: nuk")
                                if (!identical(dim(object@sk), c(K, D, D)))
                                    msg <- validMsg(msg, "Wrong dimensions: sk")
                                if (!identical(names(object@nuk), names(object@lambdak)))
                                    msg <- validMsg(msg, "nuk and lambdak names don't match")
                                if (is.null(names(object@nuk)))
                                    msg <- validMsg(msg, "Missing names")
                                if (!identical(rownames(object@mk), rownames(object@sk)))
                                    msg <- validMsg(msg, "nmk and sk rownames don't match")
                                if (!identical(rownames(object@mk), names(object@nuk)))
                                    msg <- validMsg(msg, "rownames and names don't match")
                                if (!identical(colnames(object@mk), dimnames(object@sk)[[2]]))
                                    msg <- validMsg(msg, "mk and sk[2] colnames don't match")
                                if (!identical(colnames(object@mk), dimnames(object@sk)[[3]]))
                                    msg <- validMsg(msg, "mk and sk[3] colnames don't match")
                                if (is.null(msg)) TRUE
                                else msg
                            })

setMethod("show", "ComponentParam",
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep = "")
            cat(" method:", object@method, "\n")
            cat(" Number of components:", object@K, "\n")
            cat(" Number of samples:", object@D, "\n")
            invisible(NULL)
          })


##' @title Container for a single MCMC chain results
##'
##' @aliases class:MCMCChain MCMCChain-class MCMCChain
##' @slot n `integer(1)` indicating the number of MCMC interactions.
##' @slot K `integer(1)` indicating the number of components.
##' @slot N `integer(1)` indicating the number of proteins.
##' @slot Component `matrix(N, n)` component allocation results.
##' @slot ComponentProb `matrix(N, n, K)` component allocation probabilities.
##' @slot Outlier `matrix(N, n)` outlier allocation results.
##' @slot OutlierProb `matrix(N, n, 2)` outlier allocation probabilities.
.MCMCChain <- setClass("MCMCChain",
                       slots = c(n = "integer",
                                 K = "integer",
                                 N = "integer",
                                 Component = "matrix",
                                 ComponentProb = "array",
                                 Outlier = "matrix",
                                 OutlierProb = "array",
                                 ComponentProtein = "array",
                                 ComponentParam = "ComponentParam"),
                       validity = function(object) {
                           msg <- validMsg(NULL, NULL)
                           N <- object@N
                           n <- object@n
                           K <- object@K
                           if (!identical(nrow(object@Component), N))
                               msg <- validMsg(msg, "Wrong number of proteins in component")
                           if (!identical(nrow(object@Outlier), N))
                               msg <- validMsg(msg, "Wrong number of proteins in outlier")
                           if (!identical(ncol(object@Component), n))
                               msg <- validMsg(msg, "Wrong number of iterations in component")
                           if (!identical(ncol(object@Outlier), n))
                               msg <- validMsg(msg, "Wrong number of iterations in outlier")
                           if (!identical(rownames(object@Component), rownames(object@ComponentProb)))
                               msg <- validMsg(msg, "Component rownames don't match")
                           if (!identical(rownames(object@Outlier), rownames(object@OutlierProb)))
                               msg <- validMsg(msg, "Outlier rownames don't match")
                           if (!identical(rownames(object@Outlier), rownames(object@Component)))
                               msg <- validMsg(msg, "Proteins don't match between component and outlier")
                           if (!identical(dim(object@ComponentProb)[3], K))
                               msg <- validMsg(msg, "Wrong number of components in component probability")
                           if (is.null(msg)) TRUE
                           else msg
                       })


setMethod("show", "MCMCChain",
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep = "")
            cat(" Number of components:", object@K, "\n")
            cat(" Number of proteins:", object@N, "\n")
            cat(" Number of iterations:", object@n, "\n")
            invisible(NULL)
          })


##'
##' @slot chains `list()` containing the individual full MCMC chain
##'     results. Each element must be of class `MCMCChain`.
.MCMCChains <- setClass("MCMCChains",
                        slots = c(chains = "list"),
                        validity = function(object) {
                            msg <- validMsg(NULL, NULL)
                            cls <- sapply(object@chains,
                                          function(x) inherits(x, "MCMCChain"))
                            if (!all(cls))
                                msg <- validMsg(msg, "Not all items are MCMCchains.")
                            if (is.null(msg)) TRUE
                            else msg
                        })


setMethod("length", "MCMCChains",
          function(x) length(x@chains))

setMethod("[[", "MCMCChains",
          function(x, i, j = "missing", drop = "missing") x@chains[[i]])


setMethod("[", "MCMCChains",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " chain(s) available.")
              x@chains <- x@chains[i]
              x
          })



setMethod("show", "MCMCChains",
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep = "")
            cat(" Number of chains:", length(object), "\n")
            invisible(NULL)
          })
