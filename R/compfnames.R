## helper function for compfnames to create matrix 
.calcCompNumbers <- function(flist) {
    .len <- length(flist)
    ans <- matrix(NA, nrow = .len, ncol = 3)
    colnames(ans) <- c("common", "unique1", "unique2")
    .name <- vector("character", .len)
    for (i in 1:.len) {
        .name[i] <- flist[[i]]@name
        ans[i, 1] <- length(flist[[i]]@common)
        ans[i, 2] <- length(flist[[i]]@unique1)
        ans[i, 3] <- length(flist[[i]]@unique2) 
    }
    rownames(ans) <- .name
    return(ans)
}

## helper function for compfnames to create .FeatComp S4 object
.compStrings <- function(string1, string2, all = TRUE, name = "") 
    .FeatComp(name = ifelse(all, "all", name),
              common = .shared <- intersect(string1, string2), 
              unique1 = setdiff(string1, .shared),
              unique2 = setdiff(string2, .shared),
              all = all)


setGeneric("compfnames", function(x, y, ...) standardGeneric("compfnames"))

## .FeatComp class
.FeatComp <- setClass("FeatComp", 
                    slots = list(
                        name = "character",
                        common = "character",
                        unique1 = "character",
                        unique2 = "character",
                        all = "logical"))

setValidity("FeatComp",
            function(object) {
                msg <- validMsg(NULL, NULL)
                if (any(object@unique1 %in% object@common))
                    msg <- validMsg(msg, "'unique1' features found in 'common'")
                if (any(object@unique2 %in% object@common))
                    msg <- validMsg(msg, "'unique2' features found in 'common'")
                x <- intersect(object@unique1, object@unique2)
                if (length(x) > 0)
                    stop("'unique1' and 'unique2' are not disjoint.")
                if (is.null(msg)) TRUE
                else msg
            })

setMethod("show", "FeatComp",
          function(object) {
              cat("Object of class \"", class(object), "\"", sep = "")
              cat(", '", object@name, "' features:\n", sep = "")
              cat(" Common feature:", length(object@common), "\n")
              cat(" Unique to 1:", length(object@unique1), "\n")
              cat(" Unique to 2:", length(object@unique2), "\n")
          })

setMethod("compfnames",
          c("MSnSet", "MSnSet"),
          function(x, y, 
                   fcol1 = "markers",
                   fcol2, verbose = TRUE) {     
              if (missing(fcol2)) fcol2 <- fcol1
              if (is.null(fcol1) || is.null(fcol2)) fcol1 <- fcol2 <- NULL
              
              if (!is.null(fcol1) && !is.null(fcol2)) {
                  if (!(fcol1 %in% fvarLabels(x)))
                      stop("fcol1 not in fvarLabels(x)")
                  if (!(fcol2 %in% fvarLabels(y)))
                      stop("fcol2 not in fvarLabels(y)")
              }
              
              if (!is.logical(verbose))
                  stop("value to 'verbose' is not logical")
              
              .mC <- union(unique(fData(x)[, fcol1]),
                           unique(fData(y)[, fcol2]))
              .lenmC <- length(.mC)
              ans <- vector("list", .lenmC + 1)
              
              ## for all
              ans[[1]] <- .compStrings(featureNames(x), featureNames(y))
              
              ## for marker class
              if (!is.null(fcol1)) {                  
                  for (i in 1:.lenmC) {
                      .fn1mC <- featureNames(x)[fData(x)[,fcol1] == .mC[i]]
                      .fn2mC <- featureNames(y)[fData(y)[,fcol2] == .mC[i]]
                      ans[[i+1]] <- .compStrings(.fn1mC, .fn2mC,
                                                 all = FALSE, name = .mC[i]) 
                  }
              }              
              if (verbose) print(.calcCompNumbers(ans))              
              invisible(ans)
          })
