anyUnknown <- function(x, fcol = "markers", unknown = "unknown") 
    any(fData(x)[, fcol] == unknown)

calcCompNumbers <- function(flist) {
    ## return a matrix of common, unique1, unique2
    ## for all, org1, orgi, orgn, ...
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param obj1 
##' @param obj2 
##' @param fcol1 
##' @param fcol2 
##' @return 
##' @author Laurent Gatto
compfnames <- function(obj1, obj2, fcol1 = "markers", fcol2) {
    if (missing(fcol2)) fcol2 <- fcol1

    ## if any of fcol1/2 is NULL, only compare all features
    
    ans <- list()
    ## for all features and features in each marker class
    ## number of common, unique1 and unique2

    ## 1) union of all class lables in fData(obj1)[, fcol1] and
    ## fData(obj2)[, fcol2].
    ## 2) get intersection, unique1, unique2

    msg <- calcCompNumbers(ans)
    message(msg)
    
    invisible(ans)
}
