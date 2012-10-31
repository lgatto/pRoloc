context("Machine learning tests")

data(dunkley2006)
d2 <- d1 <- dunkley2006
fvarLabels(d2)[1] <- "xx"
  
.times <- 2
.seed <- 1

test_that("knn consistency", {
  .k <- c(3, 10)
  reg1 <- knnRegularisation(d1, fcol = "markers",
                            times = .times, k = .k,
                            seed = .seed, verbose = FALSE)
  reg2 <- knnRegularisation(d2, fcol = "xx",
                            times = .times, k = .k,
                            seed = .seed, , verbose = FALSE)
  expect_equal(reg1, reg2)  
  ans1 <- knnPrediction(d1, reg1, fcol = "markers")
  ans2 <- knnPrediction(d2, reg1, fcol = "xx")
  expect_equal(ans1, ans2, check.attributes=FALSE)
})


test_that("svm consistency", {
  .cost <- 2^seq(-2, 2, 2)
  .sigma <- 10^seq(-1, 1, 1)
  reg1 <- svmRegularisation(d1, fcol = "markers",
                            cost = .cost, sigma = .sigma,
                            times = .times, 
                            seed = .seed, verbose = FALSE)
  reg2 <- svmRegularisation(d2, fcol = "xx",
                            cost = .cost, sigma = .sigma,
                            times = .times,
                            seed = .seed, verbose = FALSE)
  expect_equal(reg1, reg2)  
  ans1 <- svmPrediction(d1, reg1, fcol = "markers")
  ans2 <- svmPrediction(d2, reg1, fcol = "xx")
  expect_equal(ans1, ans2, check.attributes=FALSE)
})

test_that("nb consistency", {
  .laplace <- c(0, 5)
  reg1 <- nbRegularisation(d1, fcol = "markers",
                           laplace = .laplace,
                           times = .times, seed = .seed,
                           verbose = FALSE)
  reg2 <- nbRegularisation(d2, fcol = "xx",
                           laplace = .laplace,
                           times = .times, seed = .seed,
                           verbose = FALSE)
  expect_equal(reg1, reg2)  
  ans1 <- nbPrediction(d1, reg1, fcol = "markers")
  ans2 <- nbPrediction(d2, reg1, fcol = "xx")
  expect_equal(ans1, ans2, check.attributes=FALSE)
})

test_that("plsda consistency", {
  .ncomp <- c(3, 10)
  reg1 <- plsdaRegularisation(d1, fcol = "markers",
                              ncomp = .ncomp,
                              times = 1, seed = .seed,
                              verbose = FALSE)
  reg2 <- plsdaRegularisation(d2, fcol = "xx",
                              ncomp = .ncomp,
                              times = 1, seed = .seed,
                              verbose = FALSE)
  expect_equal(reg1, reg2)  
  ans1 <- plsdaPrediction(d1, reg1, fcol = "markers")
  ans2 <- plsdaPrediction(d2, reg1, fcol = "xx")
  expect_equal(ans1, ans2, check.attributes=FALSE)
})


test_that("nnet consistency", {
  .decay <- 10^(c(-1, -5))
  .size <- c(5, 10)
  reg1 <- nnetRegularisation(d1, fcol = "markers",
                             decay = .decay, size = .size,
                             times = .times, seed = .seed,
                             verbose = FALSE)
  reg2 <- nnetRegularisation(d2, fcol = "xx",
                             decay = .decay, size = .size,
                             times = .times, seed = .seed,
                             verbose = FALSE)
  expect_equal(reg1, reg2)
  set.seed(.seed)
  ans1 <- nnetPrediction(d1, reg1, fcol = "markers")
  set.seed(.seed)
  ans2 <- nnetPrediction(d2, reg1, fcol = "xx")
  expect_equal(ans1, ans2, check.attributes = FALSE)
})

test_that("rf consistency", {
  .mtry <- c(2, 5, 10)
  reg1 <- rfRegularisation(d1, fcol = "markers",
                           mtry = .mtry,
                           times = .times, seed = .seed,
                           verbose = FALSE)
  reg2 <- rfRegularisation(d2, fcol = "xx",
                           mtry = .mtry,
                           times = .times, seed = .seed,
                           verbose = FALSE)
  expect_equal(reg1, reg2)
  set.seed(.seed)
  ans1 <- rfPrediction(d1, reg1, fcol = "markers")
  set.seed(.seed)
  ans2 <- rfPrediction(d2, reg1, fcol = "xx")  
  expect_equal(ans1, ans2, check.attributes = FALSE)
})

test_that("ksvm consistency", {
  .cost <- 2^seq(-1, 4, 5)
  reg1 <- ksvmRegularisation(d1, fcol = "markers",
                             cost = .cost,
                             times = .times, seed = .seed,
                             verbose = FALSE)
  reg2 <- ksvmRegularisation(d2, fcol = "xx",
                             cost = .cost,
                             times = .times, seed = .seed,
                             verbose = FALSE)
  expect_equal(reg1, reg2)
  set.seed(.seed)
  ans1 <- ksvmPrediction(d1, reg1, fcol = "markers")
  set.seed(.seed)
  ans2 <- ksvmPrediction(d2, reg1, fcol = "xx")  
  expect_equal(ans1, ans2, check.attributes = FALSE)
})


