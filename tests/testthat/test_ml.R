context("Machine learning tests")

data(dunkley2006)
d2 <- d1 <- dunkley2006
i <- which(fvarLabels(d1) == "markers")
stopifnot(length(i) == 1)
fvarLabels(d2)[i] <- "xx"

.times <- 2
.seed <- 1

test_that("knn consistency", {
    .k <- c(3, 10)
    reg1 <- knnOptimisation(d1, fcol = "markers",
                            times = .times, k = .k,
                            seed = .seed, verbose = FALSE)
    reg2 <- knnOptimisation(d2, fcol = "xx",
                            times = .times, k = .k,
                            seed = .seed, verbose = FALSE)
    expect_equal(reg1, reg2)
    ans1 <- knnClassification(d1, reg1, fcol = "markers")
    ans2 <- knnClassification(d2, reg1, fcol = "xx")
    expect_equal(ans1, ans2, check.attributes=FALSE)
})


test_that("svm consistency", {
    .cost <- 2^seq(-2, 2, 2)
    .sigma <- 10^seq(-1, 1, 1)
    reg1 <- svmOptimisation(d1, fcol = "markers",
                            cost = .cost, sigma = .sigma,
                            times = .times,
                            seed = .seed, verbose = FALSE)
    reg2 <- svmOptimisation(d2, fcol = "xx",
                            cost = .cost, sigma = .sigma,
                            times = .times,
                            seed = .seed, verbose = FALSE)
    expect_equal(reg1, reg2)
    ans1 <- svmClassification(d1, reg1, fcol = "markers")
    ans2 <- svmClassification(d2, reg1, fcol = "xx")
    expect_equal(ans1, ans2, check.attributes=FALSE)
})

test_that("nb consistency", {
    .laplace <- c(0, 5)
    reg1 <- nbOptimisation(d1, fcol = "markers",
                           laplace = .laplace,
                           times = .times, seed = .seed,
                           verbose = FALSE)
    reg2 <- nbOptimisation(d2, fcol = "xx",
                           laplace = .laplace,
                           times = .times, seed = .seed,
                           verbose = FALSE)
    expect_equal(reg1, reg2)
    ans1 <- nbClassification(d1, reg1, fcol = "markers")
    ans2 <- nbClassification(d2, reg1, fcol = "xx")
    expect_equal(ans1, ans2, check.attributes=FALSE)
})

## Removed on Thu Apr  2 10:26:19 BST 2015 to reduce checkig time
## test_that("plsda consistency", {
##     .ncomp <- c(3, 10)
##     reg1 <- plsdaOptimisation(d1, fcol = "markers",
##                               ncomp = .ncomp,
##                               times = 1, seed = .seed,
##                               verbose = FALSE)
##     reg2 <- plsdaOptimisation(d2, fcol = "xx",
##                               ncomp = .ncomp,
##                               times = 1, seed = .seed,
##                               verbose = FALSE)
##     expect_equal(reg1, reg2)
##     ans1 <- plsdaClassification(d1, reg1, fcol = "markers")
##     ans2 <- plsdaClassification(d2, reg1, fcol = "xx")
##     expect_equal(ans1, ans2, check.attributes=FALSE)
## })


test_that("nnet consistency", {
    .decay <- 10^(c(-1, -5))
    .size <- c(5, 10)
    reg1 <- nnetOptimisation(d1, fcol = "markers",
                             decay = .decay, size = .size,
                             times = .times, seed = .seed,
                             verbose = FALSE)
    reg2 <- nnetOptimisation(d2, fcol = "xx",
                             decay = .decay, size = .size,
                             times = .times, seed = .seed,
                             verbose = FALSE)
    expect_equal(reg1, reg2)
    set.seed(.seed)
    ans1 <- nnetClassification(d1, reg1, fcol = "markers")
    set.seed(.seed)
    ans2 <- nnetClassification(d2, reg1, fcol = "xx")
    expect_equal(ans1, ans2, check.attributes = FALSE)
})

test_that("rf consistency", {
    .mtry <- c(2, 5, 10)
    reg1 <- rfOptimisation(d1, fcol = "markers",
                           mtry = .mtry,
                           times = .times, seed = .seed,
                           verbose = FALSE)
    reg2 <- rfOptimisation(d2, fcol = "xx",
                           mtry = .mtry,
                           times = .times, seed = .seed,
                           verbose = FALSE)
    expect_equal(reg1, reg2)
    set.seed(.seed)
    ans1 <- rfClassification(d1, reg1, fcol = "markers")
    set.seed(.seed)
    ans2 <- rfClassification(d2, reg1, fcol = "xx")
    expect_equal(ans1, ans2, check.attributes = FALSE)
})

## Removed on Thu Apr  2 10:26:19 BST 2015 to reduce checkig time
## test_that("ksvm consistency", {
##     .cost <- 2^seq(-1, 4, 5)
##     reg1 <- ksvmOptimisation(d1, fcol = "markers",
##                              cost = .cost,
##                              times = .times, seed = .seed,
##                              verbose = FALSE)
##     reg2 <- ksvmOptimisation(d2, fcol = "xx",
##                              cost = .cost,
##                              times = .times, seed = .seed,
##                              verbose = FALSE)
##     expect_equal(reg1, reg2)
##     set.seed(.seed)
##     ans1 <- ksvmClassification(d1, reg1, fcol = "markers")
##     set.seed(.seed)
##     ans2 <- ksvmClassification(d2, reg1, fcol = "xx")
##     expect_equal(ans1, ans2, check.attributes = FALSE)
## })


## test_that("deprecated consistency", {
##   .k <- c(3, 10)
##   reg1 <- knnRegularisation(d1, fcol = "markers",
##                             times = .times, k = .k,
##                             seed = .seed, verbose = FALSE)
##   reg2 <- knnRegularisation(d2, fcol = "xx",
##                             times = .times, k = .k,
##                             seed = .seed, verbose = FALSE)
##   expect_equal(reg1, reg2)
##   ans1 <- knnPrediction(d1, reg1, fcol = "markers")
##   ans2 <- knnPrediction(d2, reg1, fcol = "xx")
##   expect_equal(ans1, ans2, check.attributes=FALSE)
##   ##
##   .cost <- 2^seq(-2, 2, 2)
##   .sigma <- 10^seq(-1, 1, 1)
##   reg1 <- svmRegularisation(d1, fcol = "markers",
##                             cost = .cost, sigma = .sigma,
##                             times = .times,
##                             seed = .seed, verbose = FALSE)
##   reg2 <- svmRegularisation(d2, fcol = "xx",
##                             cost = .cost, sigma = .sigma,
##                             times = .times,
##                             seed = .seed, verbose = FALSE)
##   expect_equal(reg1, reg2)
##   ans1 <- svmPrediction(d1, reg1, fcol = "markers")
##   ans2 <- svmPrediction(d2, reg1, fcol = "xx")
##   expect_equal(ans1, ans2, check.attributes=FALSE)
##   ##
##   .ncomp <- c(3, 10)
##   reg1 <- plsdaRegularisation(d1, fcol = "markers",
##                               ncomp = .ncomp,
##                               times = 1, seed = .seed,
##                               verbose = FALSE)
##   reg2 <- plsdaRegularisation(d2, fcol = "xx",
##                               ncomp = .ncomp,
##                               times = 1, seed = .seed,
##                               verbose = FALSE)
##   expect_equal(reg1, reg2)
##   ans1 <- plsdaPrediction(d1, reg1, fcol = "markers")
##   ans2 <- plsdaPrediction(d2, reg1, fcol = "xx")
##   expect_equal(ans1, ans2, check.attributes=FALSE)
## })
