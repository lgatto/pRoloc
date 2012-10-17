.pRolocEnv <- new.env(parent=emptyenv(), hash=TRUE)

## plotting stock colors and point chars
stockcol <-  c(brewer.pal(9, "Set1"),
               "#333333", ## grey20
               "#A021EF", ## purple
               "#008A45", ## springgreen4
               "#00008A", ## blue4
               "#FF6347") ## tomato
assign("stockcol", stockcol, envir = .pRolocEnv)
getStockcol <- function() get("stockcol", envir=.pRolocEnv)

stockpch <- c(15, 17,19, 18, 23:25, 7, 9, 13, 3:4,  8)
assign("stockpch", stockpch, envir = .pRolocEnv)
getStockpch <- function() get("stockpch", envir=.pRolocEnv)


unknowncol <- "#E7E7E7" ## grey91
assign("unknowncol", unknowncol, envir = .pRolocEnv)
getUnknowncol <- function() get("unknowncol", envir=.pRolocEnv)


unknownpch <- 21
assign("unknownpch", unknownpch, envir = .pRolocEnv)
getUnknownpch <- function() get("unknownpch", envir=.pRolocEnv)

