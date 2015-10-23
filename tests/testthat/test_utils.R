context("utils")

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
  library("pRolocdata")
  data(dunkley2006)
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
    library("pRolocdata")
    data(dunkley2006)
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
    library("pRolocdata")
    data(dunkley2006)
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
