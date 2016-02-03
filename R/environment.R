.pRolocEnv <- new.env(parent=emptyenv(), hash=TRUE)

lisacols <- c("#E41A1C", "#377EB8", "#238B45", "#FF7F00", "#FFD700", "#333333", 
              "#00CED1", "#A65628", "#F781BF", "#984EA3", "#9ACD32", "#8E0152", 
              "#00008A", "#67000D", "#FDAE6B", "#66C2A5", "#276419", "#CD8C95",
              "#6A51A3", "#EEAD0E", "#0000FF", "#9ACD32", "#CD6090", "#CD5B45", 
              "#B0C4DE", "#808000", "#8B795E", "#3F007D", "#6BAED6", "#FC9272")

assign("lisacol", lisacol, envir = .pRolocEnv)

##' @rdname getStockcol
setLisacol <- function()
    assign("stockcol", lisacol, envir = .pRolocEnv)

##' @rdname getStockcol
getLisacol <- function() get("lisacol", envir=.pRolocEnv)

## plotting stock colors and point chars
stockcol <-  c(brewer.pal(9, "Set1"),
               "#333333", ## grey20
               "#A021EF", ## purple
               "#008A45", ## springgreen4
               "#00008A") ## blue4

assign("stockcol", stockcol, envir = .pRolocEnv)

stockpch <- c(19, 1, 15, 0, 17, 2, 18, 5, 7, 9, 13, 3:4,  8)
assign("stockpch", stockpch, envir = .pRolocEnv)

unknowncol <- "#E0E0E0" ## grey88
assign("unknowncol", unknowncol, envir = .pRolocEnv)

unknownpch <- 21
assign("unknownpch", unknownpch, envir = .pRolocEnv)

##' These functions allow to get/set the default colours and point
##' character that are used when plotting organelle clusters and
##' unknown features. These values are parametrised at the session
##' level. Two palettes are available: the default and Lisa's colours.
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
##' getLisacol()
##' setLisacol()
##' ## now default colours are Lisa's
##' getStockcol()
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

##' @rdname getStockcol
setStockcolGui <- function() {
    n <- length(colours())
    i <- 26
    m <- matrix(c(1:n, rep(NA, m)),
                ncol = i, nrow = i)
    ## plotting
    image(m, col = colours(),
          xaxt = "n", yaxt = "n")
    k <- seq(0, 1, length.out = i)
    kk <- expand.grid(k, k)
    kk <- kk[1:n, ]
    ## points(kk)
    ## choosing
    identifycol <- function(x, y = NULL, n = length(x), pch = 19) {
        ## from ?identify
        k <- 1
        xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
        sel <- rep(FALSE, length(x)); res <- integer(0)
        while(sum(sel) < n) {
            ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE)
            if (!length(ans)) break
            ans <- which(!sel)[ans]
            text(x[ans], y[ans], k, cex = 1.5)
            k <- k + 1
            sel[ans] <- TRUE
            res <- c(res, ans)
        }
        res
    }
    ans <- identifycol(kk)
    ans <- col2hcl(colours()[ans])
    setStockcol(ans)
}

## -------------------------------

## Annotation 

## Default params is NULL - initialised by pRoloc::setAnnotationParams()
assign("params", NULL, envir=.pRolocEnv)

lockEnvironment(.pRolocEnv,bindings=TRUE)
