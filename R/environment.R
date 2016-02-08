.pRolocEnv <- new.env(parent=emptyenv(), hash=TRUE)

stockcol <- c("#E41A1C", "#377EB8", "#238B45", "#FF7F00", "#FFD700", "#333333",
              "#00CED1", "#A65628", "#F781BF", "#984EA3", "#9ACD32", "#B0C4DE",
              "#00008A", "#8B795E", "#FDAE6B", "#66C2A5", "#276419", "#CD8C95",
              "#6A51A3", "#EEAD0E", "#0000FF", "#9ACD32", "#CD6090", "#CD5B45",
              "#8E0152", "#808000", "#67000D", "#3F007D", "#6BAED6", "#FC9272")

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

stockpch <- c(19, 1, 15, 0, 17, 2, 18, 5, 7, 9, 13, 3:4,  8)
assign("stockpch", stockpch, envir = .pRolocEnv)

unknowncol <- "#E0E0E0" ## grey88
assign("unknowncol", unknowncol, envir = .pRolocEnv)

unknownpch <- 21
assign("unknownpch", unknownpch, envir = .pRolocEnv)

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
##' getStockcol()
##' getOldcol()
getStockcol <- function() get("stockcol", envir = .pRolocEnv)

##' @param cols A vector of colour \code{characters} or \code{NULL},
##' which sets the colours to the default values.
##' @rdname getStockcol
setStockcol <- function(cols) {
    if (is.null(cols)) {
        assign("stockcol", stockcol, envir = .pRolocEnv)
    } else if (cols[1] == "lisacol") {
        setLisacol()
    } else assign("stockcol", cols, envir = .pRolocEnv)
}

##' @rdname getStockcol
getStockpch <- function() get("stockpch", envir=.pRolocEnv)

##' @param pchs A vector of \code{numeric} or \code{NULL},
##' which sets the point characters to the default values.
##' @rdname getStockcol
setStockpch <- function(pchs) {
    if (is.null(pchs)) assign("stockpch", stockpch, envir = .pRolocEnv)
    else assign("stockpch", pchs, envir = .pRolocEnv)
}

##' @rdname getStockcol
getUnknowncol <- function() get("unknowncol", envir=.pRolocEnv)


##' @param col A colour \code{character} or \code{NULL},
##' which sets the colour to \code{#E7E7E7} (\code{grey91}),
##' the default colour for unknown features.
##' @rdname getStockcol
setUnknowncol <- function(col) {
    if (is.null(col)) assign("unknowncol", unknowncol, envir = .pRolocEnv)
    else assign("unknowncol", col, envir = .pRolocEnv)
}

##' @rdname getStockcol
getUnknownpch <- function() get("unknownpch", envir=.pRolocEnv)

##' @param pch A \code{numeric} vector of length 1 or \code{NULL},
##' which sets the point character to 21, the default.
##' @rdname getStockcol
setUnknownpch <- function(pch) {
    if (is.null(pch)) assign("unknownpch", unknownpch, envir = .pRolocEnv)
    else assign("unknownpch", pch, envir = .pRolocEnv)
}

## -------------------------------

## Annotation 

## Default params is NULL - initialised by pRoloc::setAnnotationParams()
assign("params", NULL, envir=.pRolocEnv)

lockEnvironment(.pRolocEnv,bindings=TRUE)
