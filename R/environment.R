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



##' Returns the default cluster colours.
##'
##' @title Default cluster colours
##' @return A \code{character} vector.
##' @author Laurent Gatto
##' @rdname getStockcol
##' @examples
##' getStockcol()
getStockcol <- function() get("stockcol", envir=.pRolocEnv)

##' Returns the default cluster point symbols.
##'
##' @title Default cluster point symbols
##' @return A \code{numeric} vector.
##' @rdname getStockcol
##' @examples
##' getStockpch()
getStockpch <- function() get("stockpch", envir=.pRolocEnv)

##' Returns the default colour for unknown features.
##'
##' @title Default colour for unknown features
##' @return A \code{character} vector or length 1.
##' @rdname getStockcol
##' @examples
##' getUnknowncol()
getUnknowncol <- function() get("unknowncol", envir=.pRolocEnv)

##' Returns the default point character for unknown features.
##'
##' @title Default point character for unknown features
##' @return A \code{numeric} vector of length 1.
##' @rdname getStockcol
##' @examples
##' getUnknownpch()
getUnknownpch <- function() get("unknownpch", envir=.pRolocEnv)
