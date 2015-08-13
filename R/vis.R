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


SpatProtVis <- function(x, methods, methargs, ...) {
    if (missing(methods)) {
        methods <- plot2Dmethods
        methodignore <- c("scree")
        methods <- methods[!methods %in% methodignore]
    }
    if (!missing(methargs)) 
        stopifnot(length(methods) == length(methargs))
    ## TODO: manage plot2D args for the different methods, possibly as
    ## a list of metargs
    vismats <- lapply(seq_along(methods), 
                      function(i) {
                          m <- methods[i]
                          message("Producting ", m, " visualisation...")
                          args <- methargs[[i]]
                          ## FIXME
                          args <- list(object = x, plot = FALSE, args, method = m)
                          suppressMessages(do.call(plot2D, args))
                      })
    names(vismats) <- methods
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
