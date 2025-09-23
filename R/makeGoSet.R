##' @noRd
##'
##' @export
getGOFromFeatures <- function(id,
                              namespace = "cellular_component",
                              evidence = NULL,
                              params = NULL,
                              verbose = FALSE,
                              nmax = 500) {
    .Defunct()
    # if (inherits(id, "MSnSet"))
    #     id <- featureNames(id)
    # namespace <- tolower(namespace)
    # if (namespace == "all")
    #     namespace <- getGONamespace()
    # namespace <- match.arg(namespace,
    #                        getGONamespace(),
    #                        several.ok = TRUE)
    # if (is.null(params))
    #     params <- getAnnotationParams()
    # attrs <- chooseAttributesX(params)
    # if (is.null(nmax)) {
    #     ans <- getBM(attributes = c(params@filter, attrs, "go_linkage_type"),
    #                  filters = params@filter, values = id,
    #                  mart = params@mart, verbose = verbose)
    # } else {
    #     ## split filter values of there are too many
    #     ## see https://support.bioconductor.org/p/86358/
    #     from <- seq(1, length(id), by = nmax)
    #     to <- c(from[-1] - 1, length(id))
    #     ans <- vector("list", length = length(from))
    #     for (k in seq_along(from))  {
    #         .id <- id[from[k]:to[k]]
    #         ans[[k]] <- getBM(attributes = c(params@filter, attrs, "go_linkage_type"),
    #                           filters = params@filter, values = .id,
    #                           mart = params@mart, verbose = verbose)
    #     }
    #     ans <- do.call(rbind, ans)
    # }
    # sel <- ans[, attrs[2]] %in% namespace
    # if (!is.null(evidence)) {
    #     evidence <- toupper(evidence)
    #     evidence <- match.arg(evidence, getGOEvidenceCodes(),
    #                           several.ok = TRUE)
    #     if ("EXP" %in% evidence)
    #         evidence <- unique(c(evidence, "EXP", "IDA", "IPI", "IMP",
    #                              "IGI", "IEP"))
    #     if ("ISS" %in% evidence)
    #         evidence <- unique(c(evidence, "ISO", "ISA", "ISM", "IGC",
    #                              "IBA", "IBD", "IKR", "IRD", "RCA"))
    #     sel <- sel & ans$go_linkage_type %in% evidence
    # }
    # return(ans[sel,])
}


##' @noRd
##'
##' @export
makeGoSet <- function(object, params,
                      namespace = "cellular_component",
                      evidence = NULL) {
    .Defunct()
    # namespace <- match.arg(namespace,
    #                        getGONamespace(),
    #                        several.ok = TRUE)
    # if (inherits(object, "MSnSet")) {
    #     fn <- featureNames(object)
    # } else if (is.character(object)) {
    #     fn <- object
    # } else {
    #     stop("object must be an MSnSet or a character")
    # }
    # if (missing(evidence)) evidence <- NULL
    # if (missing(params))
    #     params <- getAnnotationParams()
    # if (is.null(params))
    #     stop("Please set your annotation parameters. See ?AnnotationParams for details.")
    # godf <- getGOFromFeatures(fn,
    #                           params = params,
    #                           namespace = namespace,
    #                           evidence = evidence)
    # attrs <- chooseAttributesX(params)
    # l <- lapply(fn,
    #             function(x) godf[godf[, params@filter] == x, attrs[1]])
    # allgo <- unique(unlist(l))
    # gomat <- matrix(0, length(fn), length(allgo))
    # rownames(gomat) <- fn
    # colnames(gomat) <- allgo
    # for (i in 1:nrow(gomat))
    #     gomat[i, l[[i]]] <- 1
    # goset <- new("MSnSet",
    #              exprs = gomat)
    # if (inherits(object, "MSnSet"))
    #     fData(goset) <- fData(object)
    # msg <- paste0("Constructed GO set using ", namespace, " namespace")
    # goset <- MSnbase:::logging(goset, msg)
    # if (validObject(goset))
    #     return(goset)
}
