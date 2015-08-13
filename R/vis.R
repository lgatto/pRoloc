## A class for spatial proteomics visualisation, that upon
## instantiation, pre-computes all available visualisations. Also
## stores the data annotation (spatial clusters labels)

.SpatProtVis <-
    setClass("SpatProtVis",
             slots = c(
                 vismats = "list",
                 data = "MSnSet",
                 objname = "character"),
             validity = function(object) {
                 msg <- validMsg(NULL, NULL)
                 if (!validObject(object@data))
                     msg <- validMsg(msg, "MSnSet is not valid.")                 
                 nrowval <- sapply(object@vismats,
                                   function(x) nrow(x) == nrow(object@data))
                 if (!all(nrowval))
                     msg <- validMsg(msg, "Number of features do not match.")
                 fnval <- sapply(object@vismats,
                                 function(x) identical(rownames(x),
                                                       featureNames(object@data)))
                 if (!all(fnval))
                     msg <- validMsg(msg, "Feature names do not match.")
                 if (is.null(msg)) TRUE
                 else msg
             })


SpatProtVis <- function(x, ...) {
    methodignore <- c("scree", "t-SNE")
    .plot2Dmethods <- plot2Dmethods[!plot2Dmethods %in% methodignore]
    ## TODO: manage plot2D args for the different methods, possibly as
    ## a list of metargs
    vismats <- lapply(.plot2Dmethods,
                      function(m) {
                          message("Producting ", m, " visualisation...")
                          suppressMessages(plot2D(x, method = m, plot = FALSE))
n                      })
    names(vismats) <- .plot2Dmethods
    objname <- MSnbase:::getVariableName(match.call(), "x")
    .SpatProtVis(vismats = vismats, objname = objname, data = x)
}

setMethod("show", "SpatProtVis",
          function(object) {
              cat("Object of class \"", class(object), "\"\n", sep="")
              cat(" Data:", object@objname, "\n")
              cat(" Visualisation methods:",
                  paste(names(object@vismats), collapse = ", "))
              cat("\n")
          })

setMethod("plot", c("SpatProtVis", "missing"),
          function(x, y, legend, ...) {
              for (m in names(x@vismats)) {
                  plot2D(x@data, method = x@vismats[[m]], main = m, ...)
                  if (!missing(legend))
                      addLegend(x@data, where = legend, ...)
                  if (interactive())
                      readline("Type <Return> for next plot: ")
              }
              cat("Done.\n")
          })
