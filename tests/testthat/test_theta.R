context("thetaFunctions")

## test_that("results from thetaOptimisation with all primary/auxiliary is the same as knnClassification with primary/auxiliary", {
##   library("pRolocdata")
##   library("sampling")
  
##   ## Also see function sampleMSnSet
##   data(andy2011)
##   data(andy2011goCC)
##   P <- pRoloc2:::sampleMSnSet(andy2011, seed = 1)
##   A <- pRoloc2:::sampleMSnSet(andy2011goCC, seed = 1)
  
##   relabel <- function(object, fcol = "markers", 
##                       size.validation = .2, seed) { 
##     if (!missing(seed)) {
##       seed <- as.integer(seed)
##       set.seed(seed)
##     } 
##     P <- markerSet(object, fcol)
##     data <- pRoloc:::subsetAsDataFrame(P, fcol) # NB: this function renames your fcol as "markers"
##     colnames(data)[which(colnames(data)=="markers")] <- fcol
##     ## Select validation set
##     .size <- ceiling(table(data[ ,fcol]) * size.validation)
##     .size <- .size[unique(data[ ,fcol])] 
##     validation.idxP <- strata(data, fcol, size = .size, 
##                               method = "srswor")$ID_unit  
##     validation.names <- rownames(data)[validation.idxP]
##     validation.P <- P[validation.names, ]
##     train.P <- P[-validation.idxP, ]
##     fData(train.P)$test <- as.character(fData(train.P)[, fcol])
##     fData(validation.P)$test <- rep("unknown", nrow(validation.P))
##     allP <- combine(train.P, validation.P)
##   }
    
##   P <- relabel(P, seed = 1)
##   A <- relabel(A, seed = 1)
  
##   cl <- getClasses(P, verbose=FALSE)
##   expect_equal(fData(P)$test, fData(A)$test)
  
##   ## Check for using all primary
##   theta <- rep(1, length(cl))
##   resP.th <- thetaClassification(P, A, fcol="test", theta, k=c(3,3))
##   resP.knn <- knnClassification(P, fcol="test", k=3)
  
##   expect_equal(fData(resP.th)$theta, fData(resP.knn)$knn)
##   expect_equal(fData(resP.th)$theta.scores, fData(resP.knn)$knn.scores)
  
##   ## Now test with all auxiliary data
##   theta <- rep(0, length(cl))
  
##   ## filterBinMSnSet(auxiliaryData, t=0)
  
##   resA.th <- thetaClassification(P, A, fcol="test", theta, k=c(3,3))
##   resA.knn <- knnClassification(filterBinMSnSet(A, t=0), fcol="test", k=3)
  
##   expect_equal(fData(resA.th)$theta, fData(resA.knn)$knn)
##   expect_equal(fData(resA.th)$theta.scores, fData(resA.knn)$knn.scores)
## })
