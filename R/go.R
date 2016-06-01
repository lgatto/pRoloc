##' Converts GO identifiers to/from GO terms, either explicitly or by
##' checking if (any items in) the input contains \code{"GO:"}. 
##'
##' @title Convert GO ids to/from terms
##' @param x A \code{character} of GO ids or terms.
##' @param names Should a named character be returned? Default is
##' \code{TRUE}.
##' @param keepNA Should any GO term/id names that are missing or obsolete
##' be replaced with a \code{NA}? Default is \code{TRUE}. If 
##' \code{FALSE} then the GO term/id names is kept.
##' @return A \code{character} of GO terms (ids) if \code{x} were ids
##' (terms).
##' @author Laurent Gatto
##' @examples
##' goIdToTerm("GO:0000001")
##' goIdToTerm("GO:0000001", names = FALSE)
##' goIdToTerm(c("GO:0000001", "novalid"))
##' goIdToTerm(c("GO:0000001", "GO:0000002", "notvalid"))
##' goTermToId("mitochondrion inheritance")
##' goTermToId("mitochondrion inheritance", name = FALSE)
##' goTermToId(c("mitochondrion inheritance", "notvalid"))
##' prettyGoTermId("mitochondrion inheritance")
##' prettyGoTermId("GO:0000001")
##' flipGoTermId("mitochondrion inheritance")
##' flipGoTermId("GO:0000001")
##' flipGoTermId("GO:0000001", names = FALSE)
goIdToTerm <- function(x, names = TRUE, keepNA = TRUE) {
    stopifnot(requireNamespace("GO.db"))
    stopifnot(requireNamespace("AnnotationDbi"))
    ans <- rep(NA_character_, length(x))
    names(ans) <- x
    ids <- AnnotationDbi::GOID(GO.db::GOTERM)
    i <- match(x, ids)
    k <- which(!is.na(i))
    res <- AnnotationDbi::Term(GO.db::GOTERM[i[k]])
    ans[k] <- res
    if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
    if (!names) names(ans) <- NULL
    return(ans)
}

##' @rdname goIdToTerm
goTermToId <- function(x, names = TRUE, keepNA = TRUE) {
    stopifnot(requireNamespace("GO.db"))
    stopifnot(requireNamespace("AnnotationDbi"))
    ans <- rep(NA_character_, length(x))
    names(ans) <- x
    terms <- AnnotationDbi::Term(GO.db::GOTERM)
    i <- match(x, terms)
    k <- which(!is.na(i))
    res <- AnnotationDbi::GOID(GO.db::GOTERM[i[k]])
    ans[k] <- res
    if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
    if (!names) names(ans) <- NULL
    return(ans)
}

##' @rdname goIdToTerm
flipGoTermId <- function(x, names = TRUE, keepNA = TRUE) {
    isId <- grepl("GO:", x)
    if (any(isId)) ans <- goIdToTerm(x, names, keepNA)
    else ans <- goTermToId(x, names, keepNA)
    return(ans)
}

##' @rdname goIdToTerm
prettyGoTermId <- function(x) {
    y <- flipGoTermId(x)
    if (any(grepl("GO:", x))) ans <- paste0(y, " (", x, ")")
    else ans <- paste0(x, " (", y, ")")
    return(ans)
}
