.pRolocEnv <- new.env(parent=emptyenv(), hash=TRUE)

## plotting stock colors and point chars
stockcol <-  c(brewer.pal(9, "Set1"),
               "#333333", ## grey20
               "#A021EF", ## purple
               "#008A45", ## springgreen4
               "#00008A", ## blue4
               "#FF6347") ## tomato
assign("stockcol", stockcol, envir = .pRolocEnv)

stockpch <- c(15, 17,19, 18, 23:25, 7, 9, 13, 3:4,  8)
assign("stockpch", stockpch, envir = .pRolocEnv)

unknowncol <- "#E7E7E7" ## grey91
assign("unknowncol", unknowncol, envir = .pRolocEnv)

unknownpch <- 21
assign("unknownpch", unknownpch, envir = .pRolocEnv)

##' These functions allow to get/set the default colours and
##' point character that are used when plotting organelle clusters
##' and unknown features. These values are parametrised at the
##' session level.
##'
##' @title Manage default colours and point characters
##' @return A \code{character} vector.
##' @author Laurent Gatto
##' @rdname getStockcol
##' @examples
##' ## defaults for clusters
##' getStockcol()
##' getStockpch()
##' ## unknown features
##' getUnknownpch()
##' getUnknowncol()
##' ## an example
##' library(pRolocdata)
##' data(dunkley2006)
##' par(mfrow = c(2, 1))
##' plot2D(dunkley2006, fcol = "markers", main = 'Default colours')
##' setUnknowncol("black")
##' plot2D(dunkley2006, fcol = "markers", main = 'setUnknowncol("black")')
##' getUnknowncol()
##' setUnknowncol(NULL)
##' getUnknowncol()
getStockcol <- function() get("stockcol", envir=.pRolocEnv)

##' @param cols A vector of colour \code{characters} or \code{NULL},
##' which sets the colours to the default values.
##' @return Invisibly returns \code{cols}.
##' @rdname getStockcol
setStockcol <- function(cols) {
  if (is.null(cols)) assign("stockcol", stockcol, envir = .pRolocEnv)
  else assign("stockcol", cols, envir = .pRolocEnv)
}

##' @return A \code{numeric} vector.
##' @rdname getStockcol
getStockpch <- function() get("stockpch", envir=.pRolocEnv)

##' @param pchs A vector of \code{numeric} or \code{NULL},
##' which sets the point characters to the default values.
##' @return Invisibly returns \code{pchs}.
##' @rdname getStockcol
setStockpch <- function(pchs) {
  if (is.null(pchs)) assign("stockpch", stockpch, envir = .pRolocEnv)
  else assign("stockpch", pchs, envir = .pRolocEnv)
}

##' @return A \code{character} vector or length 1.
##' @rdname getStockcol
getUnknowncol <- function() get("unknowncol", envir=.pRolocEnv)


##' @param col A colour \code{character} or \code{NULL},
##' which sets the colour to \code{#E7E7E7} (\code{grey91}),
##' the default colour for unknown features.
##' @return Invisibly returns \code{col}.
##' @rdname getStockcol
setUnknowncol <- function(col) {
  if (is.null(col)) assign("unknowncol", unknowncol, envir = .pRolocEnv)
  else assign("unknowncol", col, envir = .pRolocEnv)
}

##' @return A \code{numeric} vector of length 1.
##' @rdname getStockcol
getUnknownpch <- function() get("unknownpch", envir=.pRolocEnv)

##' @param pch A \code{numeric} vector of length 1 or \code{NULL},
##' which sets the point character to 21, the default.
##' @return Invisibly returns \code{pch}.
##' @rdname getStockcol
setUnknownpch <- function(pch) {
  if (is.null(pch)) assign("unknownpch", unknownpch, envir = .pRolocEnv)
  else assign("unknownpch", pch, envir = .pRolocEnv)
}
