context("Calculate class weights")

data(dunkley2006)

test_that("classWeights" {
    w <- 1/c(14, 45, 28, 55, 20, 46, 19, 13, 21)
    expect_true(all(classWeights(dunkley2006) == w))
})
