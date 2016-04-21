context("utils")
library("pRolocdata")
data(dunkley2006)

test_that("subsetAsDataFrame preserves col/rownames", {
    library("MSnbase")
    set.seed(1)
    data(itraqdata)
    itraqdata <- itraqdata[1:10, ]
    msnset <- quantify(itraqdata ,method = "trap",
                       reporters = iTRAQ4, verbose = FALSE)
    fData(msnset)$markers <- sample(LETTERS[1:5],
                                    nrow(itraqdata),
                                    replace=TRUE)
    .fcol <- "markers"
    dfr <- pRoloc:::subsetAsDataFrame(msnset, fcol = .fcol)
    expect_equal(featureNames(msnset), rownames(dfr))
    expect_equal(c(sampleNames(msnset), .fcol), colnames(dfr))
})


test_that("getPredictions with different thresholds", {
  ## using parameters as estimated in vignette
  ## (although these are from markers, not pd.markers)
  res <- svmClassification(dunkley2006, fcol = "pd.markers", sigma = 0.1, cost = 0.5)
  classres <- as.character(fData(res)$svm)
  allpreds <- getPredictions(res, fcol = "svm", mcol = "pd.markers", 
                             t = 0, verbose = FALSE)
  allpreds <- fData(allpreds)$svm.pred
  expect_equal(classres, allpreds)
  preds <- getPredictions(res, fcol = "svm", mcol = "pd.markers",
                          t = 1, verbose = FALSE)
  preds <- fData(preds)$svm.pred
  mrkrs <- getMarkers(res, "pd.markers", verbose = FALSE)
  names(mrkrs) <- NULL
  expect_equal(preds, mrkrs)
  tserr <- ts <- tapply(fData(res)$svm.scores, fData(res)$svm, median)
  names(tserr)[1] <- "wrongnames"  
  expect_error(getPredictions(res, fcol = "svm", mcol = "pd.markers",
                              t = tserr, verbose = FALSE))
  preds <- getPredictions(res, fcol = "svm", mcol = "pd.markers",
                          t = ts, verbose = FALSE)
  preds <- fData(preds)$svm.pred
  for (k in unique(classres)) {
    ## all thresholded predictions must be > that class score
    i <- preds == k
    expect_true(all(fData(res)[i, "svm.scores"] >= ts[k]))
    ## predicted but not in thresholded predictions must have
    ## lower scores than class scores
    j <- allpreds == k
    expect_true(all(fData(res)[j & !i, "svm.scores"] < ts[k]))
  }
})

test_that("subsetAsDataFrame keeping colnames", {
    fcol <- "new"
    cn0 <- pRoloc:::subsetAsDataFrame(dunkley2006,
                                      fcol,
                                      keepColNames = TRUE)
    cn1 <- pRoloc:::subsetAsDataFrame(dunkley2006,
                                      fcol,
                                      keepColNames = FALSE)
    cn1 <- colnames(cn1)
    cn0 <- colnames(cn0)
    n <- ncol(dunkley2006) + 1
    expect_equal(cn0[n], fcol)
    expect_equal(cn1[n], "markers")
})

test_that("fDataToUnknown function", {
    x <- getMarkers(dunkley2006, "markers", verbose = FALSE)
    ## replacing unknown by unassigned
    xx <- fDataToUnknown(dunkley2006,
                         from = "unknown", to = "unassigned")
    y <- getMarkers(xx, "markers", verbose = FALSE)
    expect_identical(x == "unknown", y == "unassigned")
    expect_identical(x != "unknown", y != "unassigned")
    ## defaults, should not affect the data
    xx <- fDataToUnknown(dunkley2006)
    y <- getMarkers(xx, "markers", verbose = FALSE)
    expect_identical(x, y)
    ## a non-existing pattern, should not affect the data
    xx <- fDataToUnknown(dunkley2006,
                         from = "qqqqqqqqqqqqqq",
                         to = "unassigned")
    y <- getMarkers(xx, "markers", verbose = FALSE)
    expect_identical(x, y)
})


test_that("fDataToUnknown function, 2nd test", {
    tmp <- dunkley2006
    fData(tmp) <-
        fData(dunkley2006)[, c("assigned", "markers.orig", "markers")]
    res1 <- fDataToUnknown(tmp, from = "ER", to = "TEST", fcol = "assigned")
    res1 <- fDataToUnknown(res1, from = "ER", to = "TEST", fcol = "markers.orig")
    res1 <- fDataToUnknown(res1, from = "ER", to = "TEST", fcol = "markers")
    res2 <- fDataToUnknown(tmp, from = "ER", to = "TEST", fcol = NULL)
    expect_identical(res1, res2)
})

test_that("anyUnknowns", {
    expect_true(pRoloc:::anyUnknown(dunkley2006))
    expect_false(pRoloc:::anyUnknown(markerMSnSet(dunkley2006)))
})

test_that("isBinary", {
    tmp <- dunkley2006
    expect_false(pRoloc:::isBinary(tmp))
    sel <- exprs(tmp) < 0.5
    exprs(tmp)[sel] <- 0
    exprs(tmp)[!sel] <- 1
    expect_true(pRoloc:::isBinary(tmp))
})

test_that("check[Sorted]FeatureNames", {
    tmp2 <- tmp <- dunkley2006
    k <- sample(nrow(tmp2))
    tmp2 <- tmp[k, ]
    expect_false(pRoloc:::checkFeatureNames(tmp, tmp2))
    expect_true(pRoloc:::checkFeatureNames(tmp, tmp2[order(k), ]))
    expect_true(pRoloc:::checkSortedFeatureNames(tmp, tmp2))
})

test_that("checkFeatureNamesOverlap", {
    o <- checkFeatureNamesOverlap(dunkley2006[1:100, ], dunkley2006[90:200, ])    
    expect_length(o$markersXY, 0)
    expect_length(o$markersX, 49)
    expect_length(o$markersY, 34)
    expect_length(o$unknownsXY, 11)
    expect_length(o$unknownsX, 40)
    expect_length(o$unknownsY, 66)
    expect_identical(sort(featureNames(markerMSnSet(dunkley2006[1:100, ]))),
                     sort(c(o$markersXY, o$markersX)))
    expect_identical(sort(featureNames(markerMSnSet(dunkley2006[90:200, ]))),
                     sort(c(o$markersXY, o$markersY)))
    expect_identical(sort(featureNames(unknownMSnSet(dunkley2006[1:100, ]))),
                     sort(c(o$unknownsXY, o$unknownsX)))
    expect_identical(sort(featureNames(unknownMSnSet(dunkley2006[90:200, ]))),
                     sort(c(o$unknownsXY, o$unknownsY)))
})

test_that("checkFvarOverlap", {
    o <- pRoloc:::checkFvarOverlap(tmp, tmp)
    expect_true(all(o$lower.mismatches == 0))
    expect_true(all(o$upper.mismatches == 0))
    tmp <- table(getMarkers(dunkley2006, verbose = FALSE))
    tmp2 <- structure(as.numeric(tmp), names = names(tmp))
    expect_equal(o$matches, tmp2)
})
