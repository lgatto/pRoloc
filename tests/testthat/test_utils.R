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
    data(tan2009r1)
    ## using parameters as estimated in vignette
    ## (although these are from markers, not pd.markers)
    res <- svmClassification(dunkley2006, fcol = "pd.markers", sigma = 0.1, cost = 0.5)
    classres <- as.character(fData(res)$svm)
    allpreds <- getPredictions(res, fcol = "svm", t = 0, verbose = FALSE)
    expect_equal(classres, allpreds)
    preds <- getPredictions(res, fcol = "svm", t = 1, verbose = FALSE)
    mrkrs <- getMarkers(res, "pd.markers", verbose = FALSE)
    names(mrkrs) <- NULL
    expect_equal(preds, mrkrs)
    tserr <- ts <- tapply(fData(res)$svm.scores, fData(res)$svm, median)
    names(tserr)[1] <- "wrongnames"
    expect_error(getPredictions(res, fcol = "svm", t = tserr, verbose = FALSE))
    preds <- getPredictions(res, fcol = "svm", t = ts, verbose = FALSE)
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


