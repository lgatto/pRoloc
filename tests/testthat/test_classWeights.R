context("Calculate class weights")

test_that("classWeights", {
    data(dunkley2006)
    w <- table(fData(dunkley2006)[, "markers"])
    w <- 1/w[names(w) != "unknown"]
    expect_identical(classWeights(dunkley2006, "markers"), w)

    data(tan2009r1)
    w <- table(fData(tan2009r1)[, "pd.markers"])
    w <- 1/w[names(w) != "unknown"]
    expect_identical(classWeights(tan2009r1, "pd.markers"), w)
})
