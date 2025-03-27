.pRolocEnv <- new.env(parent=emptyenv(), hash=TRUE)

stockcol <- c("#E41A1C", "#377EB8", "#309C17", "#FF7F00", "#FFD700", "#00CED1",
              "#A65628", "#F781BF", "#984EA3", "#9ACD32", "#B0C4DE", "#00008A",
              "#FDAE6B", "#EBB7BE", "#3F8F8F", "#CF9802", "#6A51A3", "#21E8AC",
              "#0000FF", "#1D7A3E", "#BF2A6B", "#CD5B45", "#808000", "#F21D56",
              "#67000D", "#7A0C79", "#93EDF5", "#A66A6A", "#0E438A", "#DBBCF7")

assign("stockcol", stockcol, envir = .pRolocEnv)

##' @rdname getStockcol
setLisacol <- function() message("The 'lisacol' palette is now the default option")

##' @rdname getStockcol
getLisacol <- function() message("The 'lisacol' palette is now the default option")

## plotting stock colors and point chars
oldcol <-  c(brewer.pal(9, "Set1"),
             "#333333", ## grey20
             "#A021EF", ## purple
             "#008A45", ## springgreen4
             "#00008A") ## blue4

assign("oldcol", oldcol, envir = .pRolocEnv)

##' @rdname getStockcol
getOldcol <- function() get("oldcol", envir = .pRolocEnv)

##' @rdname getStockcol
setOldcol <- function()
  assign("stockcol", oldcol, envir = .pRolocEnv)

stockpch <- c(21:25, 19, 1, 15, 0, 17, 2, 18, 5, 7, 9, 13, 3:4,  8)
# stockpch <- c(19, 1, 15, 0, 17, 2, 18, 5, 7, 9, 13, 3:4,  8)
assign("stockpch", stockpch, envir = .pRolocEnv)

unknowncol <- "#E0E0E0" ## grey88
assign("unknowncol", unknowncol, envir = .pRolocEnv)

unknownpch <- 21
assign("unknownpch", unknownpch, envir = .pRolocEnv)

stockbg <- lighten(stockcol, amount = .2)
assign("stockbg", stockbg, envir = .pRolocEnv)

unknownbg <- lighten("#E0E0E0", amount = .4)
assign("unknownbg", unknownbg, envir = .pRolocEnv)

##' These functions allow to get/set the colours and point character
##' that are used when plotting organelle clusters and unknown
##' features. These values are parametrised at the session level. Two
##' palettes are available: the default palette (previously
##' \emph{Lisa's colours}) containing 30 colours and the old
##' (original) palette, containing 13 colours.
##'
##' @title Manage default colours and point characters
##' @return The \code{set} functions set (and invisibly returns)
##' colours. The \code{get} functions returns a \code{character}
##' vector of colours. For the \code{pch} functions, \code{numeric}s
##' rather than \code{character}s.
##' @author Laurent Gatto
##' @rdname getStockcol
##' @examples
##' ## defaults for clusters
##' getStockcol()
##' getStockbg()
##' getStockpch()
##' ## unknown features
##' getUnknowncol()
##' getUnknownbg()
##' getUnknownpch()
##' ## an example
##' library(pRolocdata)
##' data(dunkley2006)
##' par(mfrow = c(2, 1))
##' plot2D(dunkley2006, fcol = "markers", main = 'Default colours')
##' setUnknowncol("black")
##' setUnknownbg("grey")
##' plot2D(dunkley2006, fcol = "markers", 
##'       main = 'setUnknowncol("black") and setUnknownbg("grey")')
##' getUnknowncol()
##' getUnknownbg()
##' setUnknowncol(NULL)
##' setUnknownbg(NULL)
##' getUnknowncol()
##' getStockcol()
##' getOldcol()
getStockcol <- function() get("stockcol", envir = .pRolocEnv)

##' @param cols A vector of colour \code{characters} or \code{NULL},
##' which sets the colours to the default values.
##' @rdname getStockcol
setStockcol <- function(cols) {
  prevcols <- getStockcol()
  if (is.null(cols)) {
    assign("stockcol", stockcol, envir = .pRolocEnv)
  } else if (cols[1] == "lisacol") {
    setLisacol()
  } else assign("stockcol", cols, envir = .pRolocEnv)
  invisible(prevcols)
}

##' @rdname getStockcol
getStockpch <- function() get("stockpch", envir=.pRolocEnv)

##' @param pchs A vector of \code{numeric} or \code{NULL},
##' which sets the point characters to the default values.
##' @rdname getStockcol
setStockpch <- function(pchs) {
  prevpch <- getStockpch()
  if (is.null(pchs)) assign("stockpch", stockpch, envir = .pRolocEnv)
  else assign("stockpch", pchs, envir = .pRolocEnv)
  invisible(prevpch)
}

##' @rdname getStockcol
getUnknowncol <- function() get("unknowncol", envir=.pRolocEnv)


##' @param col A colour \code{character} or \code{NULL},
##' which sets the colour to \code{#E7E7E7} (\code{grey91}),
##' the default colour for unknown features.
##' @rdname getStockcol
setUnknowncol <- function(col) {
  prevcol <- getUnknowncol()
  if (is.null(col)) assign("unknowncol", unknowncol, envir = .pRolocEnv)
  else assign("unknowncol", col, envir = .pRolocEnv)
  invisible(prevcol)
}

##' @rdname getStockcol
getUnknownpch <- function() get("unknownpch", envir=.pRolocEnv)

##' @param pch A \code{numeric} vector of length 1 or \code{NULL},
##' which sets the point character to 21, the default.
##' @rdname getStockcol
setUnknownpch <- function(pch) {
  prevpch <- getUnknownpch()
  if (is.null(pch)) assign("unknownpch", unknownpch, envir = .pRolocEnv)
  else assign("unknownpch", pch, envir = .pRolocEnv)
  invisible(prevpch)
}

##' @rdname getStockcol
getStockbg <- function() get("stockbg", envir = .pRolocEnv)

##' @param bg A vector of colour \code{characters} or \code{NULL},
##' which sets the background (fill) color for the open plot symbols 
##' given by pch = 21:25. to the default values.
##' @rdname getStockcol
setStockbg <- function(bg) {
  prevcols <- getStockbg()
  if (is.null(bg)) {
    assign("stockbg", stockbg, envir = .pRolocEnv)
  } else if (bg[1] == "lisacol") {
    setLisacol()
  } else assign("stockbg", bg, envir = .pRolocEnv)
  invisible(prevcols)
}

##' @rdname getStockcol
getUnknownbg <- function() get("unknownbg", envir=.pRolocEnv)

##' @param bg A colour \code{character} or \code{NULL},
##' which sets the background (fill) colour for open plot symbols
##' given by pch = 21:25 to the default colour for unknown features.
##' @rdname getStockcol
setUnknownbg <- function(bg) {
  prevcol <- getUnknownbg()
  if (is.null(bg)) assign("unknownbg", unknownbg, envir = .pRolocEnv)
  else assign("unknownbg", bg, envir = .pRolocEnv)
  invisible(prevcol)
}

## -------------------------------

## Annotation 

## Default params is NULL - initialised by pRoloc::setAnnotationParams()
assign("params", NULL, envir=.pRolocEnv)

lockEnvironment(.pRolocEnv,bindings=TRUE)
