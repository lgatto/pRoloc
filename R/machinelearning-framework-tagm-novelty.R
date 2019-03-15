##' @title Contained for full Novelty TAGM results.
##'
##' @slot pooledNoveltyChain An instance of class `NoveltyChain` after pooling
##'     all separate chains.
##' @slot noveltyChains `list()` containing the individual full Novelty TAGM
##'     chain results in an `NoveltyChains` instance. Each element must be a
##'     valid `NoveltyChain` instance.
##' @slot tagmNewclusterProb A numeric vector containing the probability that
##'     a protein belongs to a new cluster (aka phenotype).
##' @md
##' @rdname NoveltyChains
.NoveltyChains <- setClass("NoveltyChains",
                           slots = c(pooledNoveltyChain = "NoveltyChain",
                                     noveltyChains = "list",
                                     tagmNewclusterProb = "numeric"),
                           validity = function(object) {
                               msg <- validMsg(NULL, NULL)
                               cls <- sapply(object@noveltyChains,
                                             function(x) inherits(x, "NoveltyChain"))
                               if (!all(cls))
                                   msg <- validMsg(msg, "Not all items are Novlty Chains.")
                               if (is.null(msg)) TRUE
                               else msg
                           })


##' @title Container for a single Novelty TAGM mcmc chain results
##'
##' @slot psm A matrix where each entry is the posterior probability that
##'     protein i localises with protein j
##' @slot maxPear An instance of class `MaxPear` containing results from
##'     maximising the Posterior Expected Adjusted Rand index.
##' @md
##' @rdname NoveltyChains
.NoveltyChain <- setClass("NoveltyChain",
                       slots = c(psm = "matrix",
                                 mp = "MaxPear"))

##' @title Contained for results of maximising the posterior expected adjusted
##'     rand index This class wraps the `list` ouput of `mcclust::maxpear` into
##'     an object.
##'
##' @slot cl A factor containing the organelle/phenotype to which the protein
##'     has been assigned.
##' @slot value numeric(1) containing the value of PEAR.
##' @slot method character(1) describing the approach used to maximise PEAR. See
##'     `mcclust::maxpear` for details.
##' @rdname NoveltyChains
.MaxPear <- setClass("MaxPear",
                     slots = c(cl = "factor",
                               value = "numeric",
                               method = "character"))

setMethod("show", "NoveltyChains",
          function(object) {
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat("Number of chains:", length(object@noveltyChains), "\n")
              cat("PEAR: ",  object@pooledNoveltyChain@mp@value, "\n")
              cat("New phenotypes:",
                  length(unique(grep("Phenotype", object@pooledNoveltyChain@mp@cl, value = TRUE))),
                  "\n")
          })