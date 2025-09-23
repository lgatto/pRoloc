##' @noRd
##'
##' @export
goIdToTerm <- function(x, names = TRUE, keepNA = TRUE) {
    .Defunct()
    # stopifnot(requireNamespace("GO.db"))
    # stopifnot(requireNamespace("AnnotationDbi"))
    # ans <- rep(NA_character_, length(x))
    # names(ans) <- x
    # ids <- AnnotationDbi::GOID(GO.db::GOTERM)
    # i <- match(x, ids)
    # k <- which(!is.na(i))
    # res <- AnnotationDbi::Term(GO.db::GOTERM[i[k]])
    # ans[k] <- res
    # if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
    # if (!names) names(ans) <- NULL
    # return(ans)
}

##' @noRd
##'
##' @export
goTermToId <- function(x, names = TRUE, keepNA = TRUE) {
    .Defunct()
    # stopifnot(requireNamespace("GO.db"))
    # stopifnot(requireNamespace("AnnotationDbi"))
    # ans <- rep(NA_character_, length(x))
    # names(ans) <- x
    # terms <- AnnotationDbi::Term(GO.db::GOTERM)
    # i <- match(x, terms)
    # k <- which(!is.na(i))
    # res <- AnnotationDbi::GOID(GO.db::GOTERM[i[k]])
    # ans[k] <- res
    # if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
    # if (!names) names(ans) <- NULL
    # return(ans)
}

##' @noRd
##'
##' @export
flipGoTermId <- function(x, names = TRUE, keepNA = TRUE) {
    .Defunct()
    # isId <- grepl("GO:", x)
    # if (any(isId)) ans <- goIdToTerm(x, names, keepNA)
    # else ans <- goTermToId(x, names, keepNA)
    # return(ans)
}

##' @noRd
##'
##' @export
prettyGoTermId <- function(x) {
    .Defunct()
    # y <- flipGoTermId(x)
    # if (any(grepl("GO:", x))) ans <- paste0(y, " (", x, ")")
    # else ans <- paste0(x, " (", y, ")")
    # return(ans)
}
