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
    dfr <- subsetAsDataFrame(msnset, fcol = .fcol)
    expect_equal(featureNames(msnset), rownames(dfr))
    expect_equal(c(sampleNames(msnset), .fcol), colnames(dfr))
})
