context("thetaFunctions")

test_that("theta matrices passed to knntl are the same as those output and stored in hyperparameters", {
  library("pRolocdata")
  data(andy2011)
  data(andy2011goCC)

  ## Do we get the same results with and without specifying
  ## method = "Breckels"
  seed <- sample(.Machine$integer.max, 1)
  tl1 <- knntlOptimisation(andy2011, andy2011goCC,
                           fcol = "markers.orig",
                           times = 1,
                           by = 1, k = c(3, 3),
                           BPPARAM = SerialParam(),
                           method = "Breckels",
                           seed = seed)
  tl2 <- knntlOptimisation(andy2011, andy2011goCC,
                           fcol = "markers.orig",
                           times = 1,
                           by = 1, k = c(3, 3),
                           BPPARAM = SerialParam(),
                           seed = seed)
  expect_equal(tl1, tl2)

  expect_is(combineThetaRegRes(list(tl1, tl2)), "ThetaRegRes")

  par <- getParams(tl1)
  res1 <- knntlClassification(andy2011, andy2011goCC,
                              fcol = "markers.orig",
                              bestTheta = par, k = c(3,3),
                              seed = seed)
  res2 <- knntlClassification(andy2011, andy2011goCC,
                              fcol = "markers.orig",
                              bestTheta = par, k = c(3,3),
                              seed = seed)
  expect_equal(fData(res1)$knntl.scores, fData(res2)$knntl.scores)
  expect_equal(fData(res1)$knntl, fData(res2)$knntl)

  ## Test whether passing weights using the thetas function
  ## or calculating internally produces the same results
  n <- length(getMarkerClasses(andy2011,fcol = "markers.orig"))
  weights <- thetas(nclass = n, by = 1, verbose = FALSE)

  tl3 <- knntlOptimisation(andy2011, andy2011goCC,
                           fcol = "markers.orig",
                           times = 1, th = weights,
                           by = 1, k = c(3, 3),
                           BPPARAM = SerialParam(),
                           seed = seed)
  expect_equal(tl1, tl3)

  ## Test whether Wu weights are calculated correctly i.e. the same weight
  ## per class, data-specific not class-specific
  wuweights <- matrix(rep(seq(0, 1, length.out = 4), n), ncol = n)
  tl4 <-  knntlOptimisation(andy2011, andy2011goCC,
                            fcol = "markers.orig",
                            times = 1,
                            length.out = 4,
                            k = c(3, 3),
                            BPPARAM = SerialParam(),
                            method = "Wu",
                            seed = seed)
  tl5 <-  knntlOptimisation(andy2011, andy2011goCC,
                            fcol = "markers.orig",
                            times = 1,
                            length.out = 4,
                            k = c(3, 3),
                            BPPARAM = SerialParam(),
                            th = wuweights,
                            seed = seed)
  expect_equal(tl4, tl5)

  ## Test if the auxiliary is a bigger (more proteins) *and* the markers
  ## are the same we get the same results
  torm <- sample(which(fData(andy2011)$markers.orig == "unknown"), 100)
  andysmall <- andy2011[-torm, ]

  tl1 <- knntlOptimisation(andysmall, andy2011goCC,
                           fcol = "markers.orig",
                           times = 1,
                           by = 1, k = c(3, 3),
                           BPPARAM = SerialParam(),
                           seed = seed)

  tl2 <- knntlOptimisation(andy2011, andy2011goCC,
                           fcol = "markers.orig",
                           times = 1,
                           by = 1, k = c(3, 3),
                           BPPARAM = SerialParam(),
                           seed = seed)

  expect_identical(tl1, tl2)

  ## Test if the auxiliary is a bigger (more proteins) *and* the markers
  ## are different code doesn't break
  torm <- sample(which(fData(andy2011)$markers.orig != "unknown"), 100)
  andysmall <- andy2011[-torm, ]
  tl <- knntlOptimisation(andysmall, andy2011goCC,
                          fcol = "markers.orig",
                          times = 1,
                          by = 1, k = c(3, 3),
                          BPPARAM = SerialParam(),
                          seed = seed)
  expect_true(validObject(tl))

#   ## check logging and no-logging results in the same results
#   seed <- sample(.Machine$integer.max, 1)
#   tl1 <- knntlOptimisation(andy2011, andy2011goCC,
#                            fcol = "markers.orig",
#                            times = 1,
#                            by = 1, k = c(3, 3),
#                            BPPARAM = SerialParam(),
#                            seed = seed)
#
#   tl2 <- knntlOptimisation(andy2011, andy2011goCC,
#                            fcol = "markers.orig",
#                            times = 1,
#                            by = 1, k = c(3, 3),
#                            BPPARAM = SerialParam(),
#                            seed = seed,
#                            log = TRUE)
#
#   expect_identical(tl1@results, tl2@results)
#   expect_identical(tl1@hyperparameters, tl2@hyperparameters)

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
