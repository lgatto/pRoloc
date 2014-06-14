anyUnknown <- function(x, fcol = "markers", unknown = "unknown") 
    any(fData(x)[, fcol] == unknown)

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

## .FeatComp class
.FeatComp <- setClass("FeatComp", 
                    slots = list(
                        name = "character",
                        common = "character",
                        unique1 = "character",
                        unique2 = "character",
                        all = "logical"))

## helper function for compfnames to create .FeatComp S4 object
.compStrings <- function(string1, string2, all = TRUE, name = "") {
    if (all) 
        name <- "all"
    
    .shared <- intersect(string1, string2)
    .uniq1 <- setdiff(string1, .shared)
    .uniq2 <- setdiff(string2, .shared)
    
    .FeatComp(name = name,
            common = .shared,
            unique1 = .uniq1,
            unique2 = .uniq2,
            all = all)
}

##' @description \code{compfnames} allows to compare two objects of class 
##' \code{"MSnSet"} with regard to their occurrence of features in total and 
##' separated to unique marker levels on the basis of their \code{featureNames}.
##' Calling \code{compfnames} prints a 
##' matrix giving information at first glance about similarities and 
##' differences. Invisibly a list containing the individual \code{featureNames}
##' is returned.
##' @details \code{featureNames(obj1)} and \code{featureNames(obj2)} are used 
##' to compute similarities and differences in the occurence of features. 
##' \code{fcol1} and \code{fcol2} take measured variable names or \code{NULL} as 
##' input which may be of great expedience in comparative spatial proteomics.
##' @title compfnames
##' @param obj1 an object of class \code{"MSnSet"}.
##' @param obj2 an object of class \code{"MSnSet"}.
##' @param fcol1 an object of class \code{"character"} either in 
##' \code{fvarLabels(obj1)} or \code{NULL}.
##' @param fcol2 an object of class \code{"character"} either in 
##' \code{fvarLabels(obj2)} or \code{NULL}.
##' @param verbose an object of class \code{"logical"}. Default \code{TRUE}.
##' Set to \code{FALSE} to suppress printing matrix.
##' @return \code{compfnames} prints a matrix showing features common in 
##' \code{obj1} and \code{obj2} and being unique to \code{obj1} and \code{obj2}.
##' This is done in total and for each unique marker level. Invisibly it 
##' returns a list which contains feature names split accordingly to 
##' \code{"all"} and to each unique marker level giving information about 
##' features in common and unique for \code{obj1} and \code{obj2}.
##' @author Laurent Gatto <lg390@@cam.ac.uk>, Thomas Naake <tn299@@cam.ac.uk>
##' @export
compfnames <- function(obj1, obj2, 
                            fcol1 = "markers", fcol2, verbose = TRUE) {    
    
    if (!inherits(obj1, "MSnSet"))
        stop("obj1 not of class MSnSet")
    if (!inherits(obj2, "MSnSet"))
        stop("obj2 not of class MSnSet")
    
    if (missing(fcol2)) fcol2 <- fcol1
    if (is.null(fcol1) || is.null(fcol2)) fcol1 <- fcol2 <- NULL
    
    if (!is.null(fcol1) && !is.null(fcol2)) {
        if (!(fcol1 %in% fvarLabels(obj1)))
            stop("fcol1 not in fvarLabels(obj1)")
        if (!(fcol2 %in% fvarLabels(obj2)))
            stop("fcol2 not in fvarLabels(obj2)")
    }
    
    if (!is.logical(verbose))
        stop("value to 'verbose' is not logical")
    
    .mC <- union(unique(fData(obj1)[, fcol1]), unique(fData(obj2)[, fcol2]))
    .lenmC <- length(.mC)
    ans <- vector("list", .lenmC + 1)
    
    ## for all
    ans[[1]] <- .compStrings(featureNames(obj1), featureNames(obj2))   
    
    ## for marker class
    if (!is.null(fcol1)) {
    
        for (i in 1:.lenmC) {
            .fn1mC <- featureNames(obj1)[fData(obj1)[,fcol1] == .mC[i]]
            .fn2mC <- featureNames(obj2)[fData(obj2)[,fcol2] == .mC[i]]
            ans[[i+1]] <- .compStrings(.fn1mC, .fn2mC, all = FALSE, name = .mC[i]) 
        }
    } 
    
    if (verbose)
        print(.calcCompNumbers(ans))
    
    invisible(ans)
}
