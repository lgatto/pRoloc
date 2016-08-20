##' The function pulls the gene ontology (GO) terms
##' for a set of feature names. 
##'
##' @title Retrieve GO terms for feature names
##' @param id An \code{character} with feature names to be pulled from
##'     biomart. If and \code{MSnSet} is provided, then
##'     \code{featureNames(id)} is used.
##' @param namespace The GO namespace. One of
##'     \code{biological_process}, \code{cellular_component} (default)
##'     or \code{molecular_function}.
##' @param evidence The GO evidence code. See
##'     \code{showGOEvidenceCodes} for details. If \code{NULL}
##'     (default), no filtering based on the evidence code is
##'     performed.
##' @param params An instance of class
##'     \code{"\linkS4class{AnnotationParams}"}.
##' @param verbose A \code{logical} defining verbosity of the
##'     function. Default is \code{FALSE}.
##' @param nmax As described in
##'     \url{https://support.bioconductor.org/p/86358/}, the Biomart
##'     result can be unreliable for large queries. This argument
##'     splits the input in chunks of length \code{nmax} (default is
##'     1000). If set to \code{NULL}, the query is performed in full.
##' @return A \code{data.frame} with relevant GO terms.
##' @author Laurent Gatto
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' data(dunkley2006params)
##' dunkley2006params
##' fn <- featureNames(dunkley2006)[1:5]
##' getGOFromFeatures(fn, params = dunkley2006params)
getGOFromFeatures <- function(id,
                              namespace = "cellular_component",
                              evidence = NULL,
                              params = NULL,
                              verbose = FALSE,
                              nmax = 1000) {
    if (inherits(id, "MSnSet"))
        id <- featureNames(id)
    namespace <- tolower(namespace)
    if (namespace == "all")
        namespace <- getGONamespace()
    namespace <- match.arg(namespace,
                           getGONamespace(),
                           several.ok = TRUE)
    if (is.null(params))
        params <- getAnnotationParams()
    attrs <- chooseAttributesX(params)
    if (is.null(nmax)) {
        ans <- getBM(attributes = c(params@filter, attrs, "go_linkage_type"),
                     filters = params@filter, values = id,
                     mart = params@mart, verbose = verbose)
    } else {
        ## split filter values of there are too many
        ## see https://support.bioconductor.org/p/86358/
        from <- seq(1, length(id), by = nmax)
        to <- c(from[-1] - 1, length(id))
        ans <- vector("list", length = length(from))
        for (k in seq_along(from))  {
            .id <- id[from[k]:to[k]]
            ans[[k]] <- getBM(attributes = c(params@filter, attrs, "go_linkage_type"),
                              filters = params@filter, values = .id,
                              mart = params@mart, verbose = verbose)
        }
        ans <- do.call(rbind, ans)
    }
    sel <- ans[, attrs[2]] %in% namespace
    if (!is.null(evidence)) {
        evidence <- toupper(evidence)
        evidence <- match.arg(evidence, getGOEvidenceCodes(),
                              several.ok = TRUE)
        if ("EXP" %in% evidence)
            evidence <- unique(c(evidence, "EXP", "IDA", "IPI", "IMP",
                                 "IGI", "IEP"))
        if ("ISS" %in% evidence)
            evidence <- unique(c(evidence, "ISO", "ISA", "ISM", "IGC",
                                 "IBA", "IBD", "IKR", "IRD", "RCA"))
        sel <- sel & ans$go_linkage_type %in% evidence
    }
    return(ans[sel,])
}


##' Creates a new \code{"\linkS4class{MSnSet}"} instance populated
##' with a GO term binary matrix based on an original \code{object}.
##'
##' @title Creates a GO feature \code{MSnSet}
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}
##' or a character of feature names.
##' @param params An instance of class \code{"\linkS4class{AnnotationParams}"},
##' compatible with \code{featureNames(object)}'s format.
##' @param namespace The ontology name space. One or several of
##' \code{"biological_process"},  \code{"cellular_component"} or
##' \code{"molecular_function"}.
##' @param evidence GO evidence filtering.
##' @return A new \code{"\linkS4class{MSnSet}"} with the GO terms
##' for the respective features in the original \code{object}.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' data(dunkley2006params)
##' goset <- makeGoSet(dunkley2006[1:10, ],
##'                    dunkley2006params)
##' goset
##' exprs(goset)[1:10, 1:5]
##' image(goset)
makeGoSet <- function(object, params,
                      namespace = "cellular_component",
                      evidence = NULL) {
    namespace <- match.arg(namespace,
                           getGONamespace(),
                           several.ok = TRUE)
    if (inherits(object, "MSnSet")) {
        fn <- featureNames(object)
    } else if (is.character(object)) {
        fn <- object
    } else {
        stop("object must be an MSnSet or a character")
    }
    if (missing(evidence)) evidence <- NULL
    if (missing(params))
        params <- getAnnotationParams()
    if (is.null(params))
        stop("Please set your annotation parameters. See ?AnnotationParams for details.")
    godf <- getGOFromFeatures(fn,
                              params = params,
                              namespace = namespace,
                              evidence = evidence)
    attrs <- chooseAttributesX(params)
    l <- lapply(fn,
                function(x) godf[godf[, params@filter] == x, attrs[1]]) 
    allgo <- unique(unlist(l))
    gomat <- matrix(0, length(fn), length(allgo))
    rownames(gomat) <- fn
    colnames(gomat) <- allgo
    for (i in 1:nrow(gomat)) 
        gomat[i, l[[i]]] <- 1
    goset <- new("MSnSet",
                 exprs = gomat)
    if (inherits(object, "MSnSet"))
        fData(goset) <- fData(object)
    msg <- paste0("Constructed GO set using ", namespace, " namespace")
    goset <- MSnbase:::logging(goset, msg)
    if (validObject(goset))
        return(goset)  
}

