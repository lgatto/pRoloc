library("pRoloc")
library("pRolocdata")
data(dunkley2006)

## Test original PD
## dunkley2006pdunittest1 <- phenoDisco(dunkley2006, times=100, allIter=TRUE, 
##                                     seed=1, GS=10, verbose=TRUE)
## save(dunkley2006pdunittest1, file="dunkley2006pdunittest1.rda")

## Test new PD
## dunkley2006pdunittest2 <- phenoDisco2(dunkley2006, times=100, allIter=TRUE, 
##                                     seed=1, GS=10, verbose=TRUE, BPPARAM=NULL)
## save(dunkley2006pdunittest2, file="dunkley2006pdunittest2.rda")
