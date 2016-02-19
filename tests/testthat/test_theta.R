context("thetaFunctions")

test_that("theta matrices passed to knntl are the same as those output and stored in hyperparameters", {
  library("pRolocdata")
  data(andy2011)
  data(andy2011goCC)
  tl1 <- knntlOptimisation(andy2011, andy2011goCC,
                           fcol = "markers.orig",
                           times = 1,
                           by = 1, k = c(3, 3), 
                           BPPARAM = MulticoreParam(1L), 
                           method = "Breckels")
  tl2 <-knntlOptimisation(andy2011, andy2011goCC,
                          fcol = "markers.orig",
                          times = 1,
                          by = 1, k = c(3, 3), 
                          BPPARAM = MulticoreParam(1L))
  expect_equal(tl1@hyperparameters, tl2@hyperparameters)
  
  weights <- thetas(nclass = length(getMarkerClasses(andy2011, 
                                                     fcol = "markers.orig", 
                                                     verbose = FALSE)),
                    by = 1, verbose = FALSE)
  
  tl3 <- knntlOptimisation(andy2011, andy2011goCC,
                           fcol = "markers.orig",
                           times = 1, th = weights,
                           by = 1, k = c(3, 3), 
                           BPPARAM = MulticoreParam(1L))
  expect_equal(tl3@hyperparameters, tl1@hyperparameters)
  
  n <- length(getMarkerClasses(andy2011, fcol = "markers.orig", verbose = FALSE))
  wuweights <- matrix(rep(seq(0, 1, length.out = 4), n), ncol = n)
  tl4 <-  knntlOptimisation(andy2011, andy2011goCC,
                            fcol = "markers.orig",
                            times = 1,
                            length.out = 4, 
                            k = c(3, 3), 
                            BPPARAM = MulticoreParam(1L),
                            method = "Wu")
  dimnames(tl4@hyperparameters$theta) <- NULL
  expect_equal(wuweights, tl4@hyperparameters$theta)
}) 

## test_that("results from thetaOptimisation with all primary/auxiliary is the same as knnClassification with primary/auxiliary", {
##   library("pRolocdata")
##   library("sampling")
  
##   ## Also see function sampleMSnSet
##   data(andy2011)
##   data(andy2011goCC)
##   P <- sampleMSnSet(andy2011, seed = 1)
##   A <- sampleMSnSet(andy2011goCC, seed = 1)
  
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
  
##   cl <- getMarkerClasses(P, verbose=FALSE)
##   expect_equal(fData(P)$test, fData(A)$test)
  
##   ## Check for using all primary
##   theta <- rep(1, length(cl))
##   resP.th <- thetaClassification(P, A, fcol="test", theta, k=c(3,3))
##   resP.knn <- knnClassification(P, fcol="test", k=3)
  
##   expect_equal(fData(resP.th)$knntl, fData(resP.knn)$knn)
##   expect_equal(fData(resP.th)$knntl.scores, fData(resP.knn)$knn.scores)
  
##   ## Now test with all auxiliary data
##   theta <- rep(0, length(cl))
  
##   ## filterBinMSnSet(auxiliaryData, t=0)
  
##   resA.th <- thetaClassification(P, A, fcol="test", theta, k=c(3,3))
##   resA.knn <- knnClassification(filterBinMSnSet(A, t=0), fcol="test", k=3)
  
##   expect_equal(fData(resA.th)$knntl, fData(resA.knn)$knn)
##   expect_equal(fData(resA.th)$knntl.scores, fData(resA.knn)$knn.scores)
## })
