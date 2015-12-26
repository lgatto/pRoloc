context("Calculate class weights")

data(dunkley2006)

test_that("classWeights" {
    w <- 1/c(14, 45, 28, 55, 20, 46, 19, 13, 21)
    expect_true(all(classWeights(dunkley2006) == w))
    
    ## original code from the vignette
    w <- table(fData(tan2009r1)[, "pd.markers"])
    (w <- 1/w[names(w) != "unknown"])
    expect_identical(classWeights(tan2009r1, "pd.markers"), w)
})
