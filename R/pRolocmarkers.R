## pRoloc markers
## data is in inst/extdata/marker_species.rda


##' This function retrieves a list of organelle markers or, if no
##' \code{species} is provided, prints a description of available
##' marker sets. The markers can be added to and \code{MSnSet} using
##' the \code{\link{addMarkers}} function.
##'
##' The markers have been contributed by various members of the
##' Cambridge Centre for Proteomics, in particular Dr Dan Nightingale
##' for yeast, Dr Andy Christoforou and Dr Claire Mulvey for human, Dr
##' Arnoud Groen for Arabodopsis and Dr Claire Mulvey for mouse. In
##' addition, original (curated) markers from the \code{pRolocdata}
##' datasets have been extracted (see \code{pRolocdata} for details
##' and references).  Curation involved verification of publicly
##' available subcellular localisation annotation based on the
##' curators knowledge of the organelles/proteins considered and
##' tracing the original statement in the literature.
##'
##' These markers are provided as a starting point to generate
##' reliable sets of organelle markers but still need to be verified
##' against any new data in the light of the quantitative data and the
##' study conditions.
##' 
##' @title Organelle markers
##' @param species The species of interest.
##' @return Prints a description of the available marker lists if
##' \code{species} is missing or a named character with organelle
##' markers.
##' @author Laurent Gatto
##' @seealso \code{\link{addMarkers}} to add markers to an
##' \code{MSnSet} and \code{\link{markers}} for more information about
##' marker encoding.
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
            cat(" Ids: ", m$ids, ", ", length(m$markers),
                " markers\n", sep = "")
        }
    } else {
        if (species == "scer")
            species <- "scer_uniprot"
        species <- tolower(species)
        x <- sub(".rds", "", sub("^.+marker_", "", fls))
        k <- match(species, x)
        if (is.na(k))
            stop("Available species: ", paste(x, collapse = ", "),
                 ". See pRolocmarkers() for details.")
        m <- readRDS(fls[k])
        return(m$markers)
    }
}
