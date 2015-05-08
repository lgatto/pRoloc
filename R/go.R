##' Converts GO identifiers to/from GO terms, either explicitly or by
##' checking if (any items in) the input contains \code{"GO:"}. 
##'
##' @title Convert GO ids to/from terms
##' @param x A \code{character} of GO ids or terms.
##' @return A \code{character} of GO terms (ids) if \code{x} were ids
##' (terms).
##' @author Laurent Gatto
##' @examples
##' goIdToTerm("GO:0000001")
##' goIdToTerm(c("GO:0000001", "novalid"))
##' goIdToTerm(c("GO:0000001", "GO:0000002", "notvalid"))
##' goTermToId("mitochondrion inheritance")
##' goTermToId(c("mitochondrion inheritance", "notvalid"))
##' prettyGoTermId("mitochondrion inheritance")
##' prettyGoTermId("GO:0000001")
##' flipGoTermId("mitochondrion inheritance")
##' flipGoTermId("GO:0000001")
goIdToTerm <- function(x) {
    stopifnot(requireNamespace("GO.db"))
    stopifnot(requireNamespace("AnnotationDbi"))
    ans <- rep(NA_character_, length(x))
    names(ans) <- x
    ids <- AnnotationDbi::GOID(GO.db::GOTERM)
    i <- match(x, ids)
    k <- which(!is.na(i))
    res <- AnnotationDbi::Term(GO.db::GOTERM[i[k]])
    ans[k] <- res
    return(ans)
}

##' @rdname goIdToTerm
goTermToId <- function(x) {
    stopifnot(requireNamespace("GO.db"))
    stopifnot(requireNamespace("AnnotationDbi"))
    ans <- rep(NA_character_, length(x))
    names(ans) <- x
    terms <- AnnotationDbi::Term(GO.db::GOTERM)
    i <- match(x, terms)
    k <- which(!is.na(i))
    res <- AnnotationDbi::GOID(GO.db::GOTERM[i[k]])
    ans[k] <- res
    return(ans)
}

##' @rdname goIdToTerm
flipGoTermId <- function(x) {
    isId <- grepl("GO:", x)
    if (any(isId)) ans <- goIdToTerm(x)
    else ans <- goTermToId(x)
    return(ans)
}

##' @rdname goIdToTerm
prettyGoTermId <- function(x) {
    y <- flipGoTermId(x)
    if (any(grepl("GO:", x))) ans <- paste0(y, " (", x, ")")
    else ans <- paste0(x, " (", y, ")")
    return(ans)
}
