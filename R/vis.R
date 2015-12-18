setClassUnion("listOrNULL", c("list", "NULL"))

.SpatProtVis <-
    setClass("SpatProtVis",
             slots = c(
                 vismats = "list",
                 data = "MSnSet",
                 methargs = "listOrNULL",
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


SpatProtVis <- function(x, methods, dims,
                        methargs, ...) {
    if (missing(methods)) {
        methods <- plot2Dmethods
        methodignore <- c("scree")
        methods <- methods[!methods %in% methodignore]
    }
    stopifnot(methods %in% plot2Dmethods)
    if (missing(dims)) 
        dims <- replicate(length(methods), 1:2, simplify=FALSE)
    stopifnot(length(methods) == length(dims))
    if (!missing(methargs)) stopifnot(length(methods) == length(methargs))
    else methargs <- NULL
    vismats <- lapply(seq_along(methods), 
                      function(i) {
                          m <- methods[i]
                          message("Producting ", m, " visualisation...")
                          .dims <- dims[[i]]
                          if (m == "MDS" | is.null(methargs)) { 
                              suppressMessages(plot2D(x, plot = FALSE,
                                                      dims = .dims,
                                                      method = m))
                          } else {
                              .args <- methargs[i]
                              suppressMessages(plot2D(x, plot = FALSE,
                                                      dims = .dims,
                                                      method = m,
                                                      methargs = .args))
                          }
                      })
    names(vismats) <- methods
    objname <- MSnbase:::getVariableName(match.call(), "x")
    .SpatProtVis(vismats = vismats,
                 objname = objname,
                 methargs = methargs,
                 data = x)
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
              for (i in seq_along(x@vismats)) {
                  plot2D(object = x@vismats[[i]],
                         method = "none",
                         methargs = list(x@data),
                         main = names(x@vismats)[i], ...)
                  if (!missing(legend))
                      addLegend(x@data, where = legend, ...)
                  if (interactive())
                      readline("Type <Return> for next plot: ")
              }
              cat("Done.\n")
          })
