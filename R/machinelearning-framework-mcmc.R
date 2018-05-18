.MCMCChains <- setClass("MCMCChains",
                        slots = c(chains = "list"),
                        validity = function(object) {
                            msg <- validMsg(NULL, NULL)
                            cls <- sapply(object@chains,
                                          function(x) inherits(x, "MCMCchain"))
                            if (!all(cls))
                                msg <- validMsg(msg, "Not all items are MCMCchains.")
                            if (is.null(msg)) TRUE
                            else msg
                        })

setMethod("length", "MCMCChains",
          function(x) length(x@chains))

.MCMCChain <- setClass("MCMCChain",
                       slots = c(n = "integer",
                                 proteins = "list",
                                 component = "list",
                                 component_protein = "list"),
                       validity = function(object) {
                           n <- object@n
                           msg <- validMsg(NULL, NULL)
                           if (!identical(length(object@proteins), n))
                               msg <- validMsg(msg, "Protein length is not valid.")
                           if (!identical(length(object@component), n))
                               msg <- validMsg(msg, "Protein compondent is not valid.")
                           if (!identical(length(object@component_protein), n))
                               msg <- validMsg(msg, "Compoment/protein length is not valid.")
                            if (is.null(msg)) TRUE
                            else msg
                       })

.MCMCSummary <- setClass("MCMCSummary",
                         slots = c(summary = "list"))

.MCMCParams <- setClass("MCMCParams",
                        slots = c(metadata = "list",
                                  algorithm = "character",
                                  chains = "MCMCChains",
                                  summary = "MCMCSummary"))
setMethod("show", "MCMCParams",
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep = "")
            cat("Algorithm:", object@algorithm, "\n")
            cat("Number of chains:", length(object@chains), "\n")
            invisible(NULL)
          })
