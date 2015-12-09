checkAttr0 <- function(x, attr0 = pRoloc:::getAttributesOfInterest0()) 
    attr0 %in% x$name


checkAttrX <- function(x, attrX = pRoloc:::getAttributesOfInterestX()) 
    sapply(attrX,
           function(.attrX) {
               stopifnot(length(.attrX) == 2)
               .attrX[1] %in% x$name | .attrX[2] %in% x$name
           })

.MartInstance <- setClass("MartInstance",
                    slots = c(name = "character",
                              host = "character",
                              path = "character",
                              datasets = "data.frame",
                              created = "character"))

MartInstance <- function(name, host, path, mart) {
    if (missing(mart)) {
        if (missing(name)) {
            .marts <- biomaRt::listMarts(mart = NULL, host, path = path)
            i <- menu(sapply(.marts, paste, collapse = " - "))
            name <- .marts[i, "biomart"]
        }
        mart <- useMart(name, host = host, path = path)
    }
    datasets <- biomaRt::listDatasets(mart)
    datasets$dataset <- as.character(datasets$dataset)
    .MartInstance(name = name,
                  host = host,
                  path = path,
                  datasets = datasets,
                  created = date())
}

nDatasets <- function(x) nrow(x@datasets)

as.data.frame.MartInstance <- function(x) {
    ans <- x@datasets
    ans$MartInterface <- x@name
    ans
}

.MartInstanceList <- setClass("MartInstanceList",
                              slots = c(xData = "list"))

MartInstanceList <- function(x) .MartInstanceList(xData = x)

setMethod("[[", "MartInstanceList",
          function(x, i) {
              x@xData[[i]]
          })

setMethod("[", "MartInstanceList",
          function(x, i) {
              x@xData <- x@xData[i]
          })

setMethod("lapply", "MartInstanceList",
          function(X, FUN, ...) {
              ans <- lapply(X@xData, FUN, ...)
              X@xData <- ans
              X
          })

setMethod("sapply", "MartInstanceList",
          function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) 
              sapply(X@xData, FUN, ..., simplify = TRUE, USE.NAMES = TRUE))


as.data.frame.MartInstanceList <- function(x) {
    ans <- do.call("rbind", lapply(x@xData, as.data.frame.MartInstance))
    rownames(ans) <- NULL
    ans
}

setMethod("show", "MartInstance",
          function(object) {
              cat("Name:", object@name, "\n")
              cat("Created:", object@created, "\n")
              n <- nDatasets(object)
              k <- min(n, 3)
              cat("Datasets (", n, "):\n  ",
                  paste(head(object@datasets[1:k, 1]), collapse = ", "),
                  "\n", sep = "")
          })


## filter out datasets that don't have all required attributes
filterAttrs <- 
    function(x) {
        ds <- x@datasets
        n0 <- nrow(ds)
        if (n0 == 0) return(x)
        pb <- txtProgressBar(min = 0, max = nrow(ds), style = 3)
        sel <- logical(nrow(ds))
        for (i in 1:n0) {
            setTxtProgressBar(pb, i)
            .mart <- useMart(x@name,
                             dataset = ds[i, "dataset"],
                             host = x@host,
                             path = x@path)
            .attrs <- listAttributes(.mart)
            sel[i] <- all(checkAttr0(.attrs)) & all(checkAttrX(.attrs))
        }
        close(pb)
        x@datasets <- x@datasets[sel, ]
        n1 <- nrow(x@datasets)
        message(x@name, ": ", n0, " -> ", n1)
        x
    }

getMartInstanceList <- function() {
    f <- dir(system.file("extdata/", package = "pRoloc"),
             pattern = "mil.rd", full.names = TRUE)
    readRDS(file = f)
}

getMartTab <- function()     
    as.data.frame(getMartInstanceList())


getFilterList <- function(ds) {
    marttab <- getMartTab()
    ds2 <- marttab[marttab$dataset == ds, ]
    mil <- getMartInstanceList()
    i <- which(sapply(mil, slot, name = "name") == ds2[, "MartInterface"])
    mi <- mil[[i]]
    mart <- useMart(ds2[, "MartInterface"],
                    dataset = ds,
                    host = mi@host, path = mi@path)
    biomaRt::listFilters(mart)
}
