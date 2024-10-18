## pRoloc markers
## data is in inst/extdata/marker_species.rda

##' This function retrieves a list of organelle markers or, if no \code{species}
##' is provided, prints a description of available marker sets. The markers can
##' be added to and \code{MSnSet} using the \code{\link{addMarkers}}
##' function. Several marker version are provided (see Details for additional
##' information).
##'
##' @details
##'
##' Version 1 of the markers have been contributed by various members of the
##' Cambridge Centre for Proteomics, in particular Dr Dan Nightingale for yeast,
##' Dr Andy Christoforou and Dr Claire Mulvey for human, Dr Arnoud Groen for
##' Arabodopsis and Dr Claire Mulvey for mouse. In addition, original (curated)
##' markers from the \code{pRolocdata} datasets have been extracted (see
##' \code{pRolocdata} for details and references).  Curation involved
##' verification of publicly available subcellular localisation annotation based
##' on the curators knowledge of the organelles/proteins considered and tracing
##' the original statement in the literature.
##'
##' Version 2 of the markers (current default) have been updated by Charlotte
##' Hutchings from the Cambridge Centre for Proteomics. Reference species marker
##' sets are the same as those in version 1 with minor corrections and an
##' updated naming system. Version 2 also contains additional marker sets from
##' spatial proteomics publications. References for the source publications are
##' provided below:
##'
##' \itemize{
##'
##' \item Geladaki, A., Britovsek, N.K., Breckels, L.M., Smith, T.S., Vennard,
##'    O.L., Mulvey, C.M., Crook, O.M., Gatto, L. and Lilley, K.S. (2019)
##'    Combining LOPIT with differential ultracentrifugation for high-resolution
##'    spatial proteomics.  Nature Communications. 10
##'    (1). doi:10.1038/s41467-018-08191-w
##'
##' \item Christopher, J.A., Breckels, L.M., Crook, O.M., Vazquez--Chantada, M.,
##'   Barratt, D. and Lilley, K.S. (2024) Global proteomics indicates
##'   subcellular-specific anti-ferroptotic responses to ionizing
##'   radiation.p.2024.09.12.611851. doi:10.1101/2024.09.12.611851
##'
##' \item Itzhak, D.N., Tyanova, S., Cox, J. and Borner, G.H. (2016) Global,
##'   quantitative and dynamic mapping of protein subcellular localization.
##'   eLife. 5. doi:10.7554/elife.16950
##'
##' \item Villanueva, E., Smith, T., Pizzinga, M., Elzek, M., Queiroz, R.M.L.,
##'   Harvey, R.F., Breckels, L.M., Crook, O.M., Monti, M., Dezi, V., Willis,
##'   A.E. and Lilley, K.S. (2023) System-wide analysis of RNA and protein
##'   subcellular localization dynamics. Nature Methods. 1-12.
##'   doi:10.1038/s41592-023-02101-9
##'
##' \item Christoforou, A., Mulvey, C.M., Breckels, L.M., Geladaki, A., Hurrell,
##'   T., Hayward, P.C., Naake, T., Gatto, L., Viner, R., Arias, A.M. and Lilley,
##'   K.S. (2016) A draft map of the mouse pluripotent stem cell spatial
##'   proteome. Nature Communications. 7 (1). doi:10.1038/ncomms9992
##'
##' \item Barylyuk, K., Koreny, L., Ke, H., Butterworth, S., Crook, O.M.,
##'   Lassadi, I., Gupta, V., Tromer, E., Mourier, T., Stevens, T.J., Breckels,
##'   L.M., Pain, A., Lilley, K.S. and Waller, R.F. (2020) A Comprehensive
##'   Subcellular Atlas of the Toxoplasma Proteome via hyperLOPIT Provides
##'   Spatial Context for Protein Functions. Cell Host and Microbe. 28 (5),
##'   752-766.e9. doi:10.1016/j.chom.2020.09.011
##'
##' \item Moloney, N.M., Barylyuk, K., Tromer, E., Crook, O.M., Breckels, L.M.,
##'   Lilley, K.S., Waller, R.F. and MacGregor, P. (2023) Mapping diversity in
##'   African trypanosomes using high resolution spatial proteomics. Nature
##'   Communications. 14 (1), 4401. doi:10.1038/s41467-023-40125-z
##'
##' }
##'
##' Note: These markers are provided as a starting point to generate reliable
##' sets of organelle markers but still need to be verified against any new data
##' in the light of the quantitative data and the study conditions.
##'
##' @name pRolocmarkers
##'
##' @title Organelle markers
##'
##' @param species \code{character(1)} defining the species of interest. For
##'     reference species markers, this is just the species
##'     e.g. \code{"hsap"}. For published marker sets this is the species and
##'     author name e.g. \code{"hsap_geladaki"}.
##'
##' @param version \code{character(1)} defining the marker version. Default is
##'     "2".
##'
##' @return Prints a description of the available marker lists if \code{species}
##'     is missing or a named character with organelle markers.
##'
##' @author Laurent Gatto
##'
##' @seealso \code{\link{addMarkers}} to add markers to an \code{MSnSet} and
##'     \code{\link{markers}} for more information about marker encoding.
##'
##' @examples
##' pRolocmarkers()
##' pRolocmarkers("hsap")
##' table(pRolocmarkers("hsap"))
##'
##' ## Old markers
##' pRolocmarkers("hsap", version = "2")["Q9BPW9"]
##' pRolocmarkers("hsap", version = "1")["Q9BPW9"]
pRolocmarkers <- function(species, version = "2") {
    ## To add new markers:
    ##
    ## 1. Add one or multiple rds files (one per species) in inst/extdata,
    ##    making sure they all share the same prefix (say markers2_).
    ##
    ## 2. Add a new version in the supported_version character vector and the
    ##    version and prefix to the switch() statement.
    ##
    ## 3. Unless the new markers are considered experimental and not yet ready
    ##    for consumption, update the default version argument (using the
    ##    version defined in 1. above.)
    ##
    ## 4. Update the documentation to describe the new markers and provide
    ##    references.
    supported_versions <- c("1", "2")
    stopifnot(version %in% supported_versions)
    prefix <- switch(as.character(version),
                     "1" = "marker_",
                     "2" = "marker2_")
    fls <- dir(system.file("extdata", package = "pRoloc"),
               full.names = TRUE, pattern = prefix)
    if (missing(species)) {
        cat(length(fls), " marker lists (version ", version, ") available:\n",
            sep = "")
        for (f in fls) {
            m <- readRDS(f)
            x <- sub(".rds", "", sub(prefix, "", basename(f)))
            cat(m$species, " [", x, "]:\n", sep = "")
            cat(" Ids: ", m$ids, ", ", length(m$markers),
                " markers\n", sep = "")
        }
    } else {
        if (species == "scer")
            species <- "scer_uniprot"
        species <- tolower(species)
        x <- sub(".rds", "", sub(prefix, "", basename(fls)))
        k <- match(species, x)
        if (is.na(k))
            stop("Available species: ", paste(x, collapse = ", "),
                 ". See pRolocmarkers() for details.")
        m <- readRDS(fls[k])
        return(m$markers)
    }
}