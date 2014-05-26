library("pRoloc")
library("pRolocdata")
data(dunkley2006)

## Test original PD
## dunkley2006pdunittest1 <- phenoDisco(dunkley2006, times=100, allIter=TRUE, 
##                                     seed=1, GS=10, verbose=TRUE)
## save(dunkley2006pdunittest1, file="../extdata/dunkley2006pdunittest1.rda")

## Test new PD
## dunkley2006pdunittest2 <- phenoDisco2(dunkley2006, times=100, allIter=TRUE, 
##                                     seed=1, GS=10, verbose=TRUE, BPPARAM=NULL)
## save(dunkley2006pdunittest2, file="../extdata/dunkley2006pdunittest2.rda")

## Tue Apr  8 19:50:51 2014 - running pd with mclust v4.3
##  The results are different than with mclust 4.2, but arguably better
##
## dunkley2006pdunittest_20140408 <- phenoDisco(dunkley2006, times=100, allIter=TRUE, 
##                                              seed=1, GS=10, verbose=TRUE, BPPARAM=NULL)
## save(dunkley2006pdunittest_20140408, file = "../extdata/dunkley2006pdunittest_20140408.rda")
## load("../extdata/dunkley2006pdunittest2.rda")
## plot2D(dunkley2006pdunittest_20140408, fcol = "pd", fpch = "markers")
## plot2D(dunkley2006pdunittest2, fcol = "pd", fpch = "markers")

## Testing PD updates
dunkley2006pdunittest <- phenoDisco(dunkley2006, times=100,
                                    allIter=TRUE, seed=1,
                                    GS=10, verbose=TRUE, BPPARAM=NULL)

## load("../extdata/dunkley2006pdunittest2.rda")
load("../extdata/dunkley2006pdunittest_20140408.rda")
all.equal(dunkley2006pdunittest, dunkley2006pdunittest_20140408)
stopifnot(all.equal(fData(dunkley2006pdunittest), fData(dunkley2006pdunittest_20140408)))


