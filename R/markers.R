## pRoloc markers
## data is in inst/extdata/marker_species.rda


##' This function retrieves a list of organelle markers or,
##' if not \code{species} is provided, prints a description
##' of available markers sets. The markers can be added to
##' and \code{MSnSet} using the \code{\link{addMarkers}}
##' function.
##'
##' The markers have been contributed by various members of the
##' Cambridge Centre for proteomics, in particular Dan Nightingale
##' for yeast, Dr Andy Christoforou for human, Dr Arnoud Groen for
##' Arabodopsis. In addition, markers from the \code{pRolocdata} datasets
##' have been extracted. See \code{pRolocdata} for details and references.
##'
##' @title Organelle markers
##' @param species The species of interest.
##' @return Prints a description of the available marker lists if
##' \code{species} is missing or a named character with organelle
##' markers.
##' @author Laurent Gatto
##' @examples
##' pRolocmarkers()
##' table(pRolocmarkers("atha"))
##' table(pRolocmarkers("hsap"))
pRolocmarkers <- function(species) {    
    fls <- dir(system.file("extdata", package = "pRoloc"),
               full.names = TRUE, pattern = "marker_")    
    if (missing(species)) {
        cat(length(fls), "marker lists available:\n")
        for (f in fls) {
            m <- readRDS(f)
            x <- sub(".rds", "", sub("^.+marker_", "", f))        
            cat(m$species, " [", x, "]:\n", sep = "")
            cat(" Ids: ", m$ids, ", ", length(m$markers), " markers\n", sep = "")
        }
    } else {
        species <- tolower(species)
        x <- sub(".rds", "", sub("^.+marker_", "", fls))
        k <- match(species, x)
        if (is.na(k))
            stop("Available species: ", paste(x, collapse = ", "), ". See pRolocmarkers() for details.")
        m <- readRDS(fls[k])
        return(m$markers)
    }
}
