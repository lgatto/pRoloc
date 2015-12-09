setClass("AnnotationParams",
         representation = representation(
           mart = "Mart",
           martname = "character", ## can probably be removed
           dataset = "character",  ## can probably be removed
           filter = "character",
           date = "character",
           biomaRtVersion = "character"),
         ## add support for org.XX.XXX.db annotation
         contains = c("Versioned"),
         prototype = prototype(
           new("Versioned", versions = c(AnnotationParams="0.1.0"))),
         validity = function(object) {
           msg <- validMsg(NULL, NULL)
           if (!validObject(object@mart))
             msg <- validMsg("Mart object not valid.")           
           if (object@mart@biomart != object@martname)
             msg <- validMsg("Mart biomart slot and martname do not match.")
           if (object@mart@dataset != object@dataset)
             msg <- validMsg("Mart dataset slot and dataset do not match.")
           if (!(object@filter %in% object@mart@filters$name))
             msg <- validMsg("Filter not found in biomart filters.")
           if (is.null(msg)) TRUE
           else msg
         })

setMethod("show",
          signature(object="AnnotationParams"),
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep="")
            cat(" Using the '", object@martname,"' BioMart database\n", sep="")
            cat(" Using the '", object@dataset, "' dataset\n", sep="")
            cat(" Using '", object@filter, "' as filter\n", sep="")
            cat(" Created on ",object@date, "\n", sep="")
            invisible(NULL) 
          })


setAnnotationParams <- function(params = NULL,
                                inputs,
                                graphics = getOption("menu.graphics")) {
  ## TODO: allow to directly pass a Mart instance or
  ##       martname and dataset or
  ##       filter or
  ##       filter and (mart of martname+dataset)
  ##       and call the respective menus
  ## -----------------------------------------------
    if (!is.null(params)) {
        if (!inherits(params, "AnnotationParams"))
            stop("'params' must be of class 'AnnotationParams'.")
    } else {
        if (!missing(inputs)) {
            if (length(inputs) != 2)
                stop("'inputs' must contain a species and a feature type.")
            m <- getMartTab()
            spidx <- grep(inputs[1], m[, "description"])
            if (length(spidx) != 1)
                stop("Couldn't find a unique species match for '",
                     inputs[1], "'.")
            else message("Using species ", m[spidx, "description"])
            foi <- getFiltersOfInterest()
            spfilters <- getFilterList(m[spidx, "dataset"])
            spfilters <- spfilters[spfilters$name %in% foi, ]
            fltidx <- grep(inputs[2], spfilters$description)
            if (length(fltidx) != 1)
                stop("Couldn't find a unique feature type match for '",
                     inputs[2], "'.")
            else message("Using feature type ",
                         spfilters[fltidx, "description"])                        
            miname <- m[spidx, "MartInterface"]
            mil <- getMartInstanceList()
            i <- which(sapply(mil, slot, name = "name") == miname)
            mi <- mil[[i]]            
            dataset <- as.character(m[spidx, "dataset"])
            filter <- as.character(spfilters[fltidx, "name"])            
            cat("Connecting to Biomart...\n")
            mart <- useMart(miname, dataset = dataset,
                            host = mi@host,
                            path = mi@path)
            params <- new("AnnotationParams",
                          mart = mart,
                          martname = miname,
                          dataset = dataset,
                          filter = filter,
                          date = date(),
                          biomaRtVersion =
                              packageDescription("biomaRt")$Version)
        } else { ## interactiven
            if (!interactive())
                stop("Use interactive mode of provide arguments.")
            m <- getMartTab()
            sortedsp <- sort(m$description)
            spsel <- menu(sortedsp, graphics, "Select species")
            spidx <- match(sortedsp[spsel],m$description)
            
            foi <- getFiltersOfInterest()            
            ds <- m[spidx, "dataset"]
            spfilters <- getFilterList(ds)
            sortedspfilters <-
                sort(spfilters$description[spfilters$name %in% foi])
            if (length(sortedspfilters)==0)
                stop("No filters available for this species.\n",
                     "Please consider emailing the author to file a bug.")
            fltsel <- menu(sortedspfilters, graphics, "Select filter")
            fltidx <- match(sortedspfilters[fltsel],
                            spfilters$description)

            ds2 <- m[m$dataset == ds, ]
            martname <- ds2[, "MartInterface"]
            mil <- getMartInstanceList()
            i <- which(sapply(mil, slot, name = "name") == martname)
            mi <- mil[[i]]
            filter <- as.character(spfilters[fltidx, "name"])
            cat("Connecting to Biomart...\n")
            mart <- useMart(martname,
                            dataset = ds,
                            host = mi@host, path = mi@path)
            params <- new("AnnotationParams",
                          mart = mart,
                          martname = martname,
                          dataset = ds,
                          filter = filter,
                          date = date(),
                          biomaRtVersion =
                              packageDescription("biomaRt")$Version)
        }
    }
    if (validObject(params)) {
        unlockBinding("params", .pRolocEnv)
        assign("params", params, envir=.pRolocEnv)
        lockBinding("params", .pRolocEnv)
    }
    invisible(params)
}

getAnnotationParams <-
    function() get("params", envir=.pRolocEnv)


##' This function prints a textual description
##' of the Gene Ontology evidence codes.
##'
##' @title GO Evidence Codes
##' @return These functions are used for their side effects of printing
##' evidence codes and their description.
##' @author Laurent Gatto
##' @examples
##' showGOEvidenceCodes()
##' getGOEvidenceCodes()
showGOEvidenceCodes <- function() {
  cat("GO Term Evidence Code\n")
  cat(" Experimental Evidence Codes\n")
  cat("  EXP: Inferred from Experiment\n")
  cat("   IDA: Inferred from Direct Assay\n")
  cat("   IPI: Inferred from Physical Interaction\n")
  cat("   IMP: Inferred from Mutant Phenotype\n")
  cat("   IGI: Inferred from Genetic Interaction\n")
  cat("   IEP: Inferred from Expression Pattern\n")
  cat(" Computational Analysis Evidence Codes\n")
  cat("  ISS: Inferred from Sequence or Structural Similarity\n")
  cat("   ISO: Inferred from Sequence Orthology\n")
  cat("   ISA: Inferred from Sequence Alignment\n")
  cat("   ISM: Inferred from Sequence Model\n")
  cat("   IGC: Inferred from Genomic Context\n")
  cat("   IBA: Inferred from Biological aspect of Ancestor\n")
  cat("   IBD: Inferred from Biological aspect of Descendant\n")
  cat("   IKR: Inferred from Key Residues\n")
  cat("   IRD: Inferred from Rapid Divergence\n")
  cat("   RCA: inferred from Reviewed Computational Analysis\n")
  cat(" Author Statement Evidence Codes\n")
  cat("   TAS: Traceable Author Statement\n")
  cat("   NAS: Non-traceable Author Statement\n")
  cat(" Curator Statement Evidence Codes\n")
  cat("   IC: Inferred by Curator\n")
  cat("   ND: No biological Data available\n")
  cat(" Automatically-assigned Evidence Codes\n")
  cat("   IEA: Inferred from Electronic Annotation\n")
  cat(" Obsolete Evidence Codes\n")
  cat("   NR: Not Recorded\n")
}

##' @rdname showGOEvidenceCodes
getGOEvidenceCodes <- function()
    c("EXP", ## experimental evidence
      "IDA", "IPI", "IMP", "IGI", "IEP", 
      "ISS", ## computational analysis evidence
      "ISO", "ISA", "ISM", "IGC", "IBA", 
      "IBD", "IKR", "IRD", "RCA",
      ## author statement - traceable or non-traceable
      "TAS", "NAS",
      ## curator statement - inferred by of no data
      "IC", "ND",
      ##  Automatically-assigned Evidence Codes
      "IEA",
      ## obsolete - not recorded
      "NR")
