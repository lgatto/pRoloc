
##' @title Contained for full Novelty Tagm results.
##'
##' @slot pooledRes An instance of class `NoveltyChain` after pooling all separate chains.
##' @slot NoveltyChains `list()` containing the individual full Novelty Tagm chain
##'     results in an `NoveltyChains` instance. Each element must be a
##'     valid `NoveltyChain` instance.
##' @slot tagm.newcluster.prob A numeric vector containing the probability that a protein
##' belongs to a new cluster (aka phenotype).    
##' @md
##' @rdname NoveltyParams
.NoveltyChains <- setClass("NoveltyChains",
                        slots = c(pooledRes = "NoveltyChain",
                                     NoveltyChains = "list",
                                     tagm.newcluster.prob = "numeric"),
                        validity = function(object) {
                          msg <- validMsg(NULL, NULL)
                          cls <- sapply(object@NoveltyChains,
                                        function(x) inherits(x, "NoveltyChain"))
                          if (!all(cls))
                            msg <- validMsg(msg, "Not all items are Novlty Chains.")
                          if (is.null(msg)) TRUE
                          else msg
                        })


##' @title Container for a single Novelty Tagm mcmc chain results
##'
##' @slot psm A matrix where each entry is the posterior probability that
##' protein i localises with protein j
##' @slot maxPear An instance of class maxPear containing results from maximising
##' the posterior expected adjusted Rand index.
##' @md
##' @rdname NoveltyParams
.NoveltyChain <- setClass("NoveltyChain",
                       slots = c(psm = "matrix",
                                 mp = "maxPear")
                       )
##' @title Contained for results of maximising the posterior expeced adjusted rand index
##' 
##' @slot  cl The organelle/phenotype the protein has been assigned to.
##' @slot  value The value of PEaR
##' @slot  method The approach used to maximise PEAR
##' @method 
##' @rdname 
.maxPear <- setClass("maxPear",
                     slots = c(cl = "factor",
                               value = "numeric",
                               method = "character")
                     )
