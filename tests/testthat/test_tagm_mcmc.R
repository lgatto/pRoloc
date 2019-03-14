context("Tagm MCMC tests")

data(dunkley2006)
d2 <- d1 <- dunkley2006
i <- which(fvarLabels(d1) == "markers")
stopifnot(length(i) == 1)
fvarLabels(d2)[i] <- "xx"

.times <- 2
.seed <- 1

test_that("TAGM MCMC consistency", {
  .numIter <- 5
  set.seed(1)
  mcmc1 <- tagmMcmcTrain(object = d1, fcol = "markers",
                         numIter = .numIter,burnin = 0, thin = 1, numChains = 2,
                         BPPARAM = SerialParam())
  set.seed(1)
  mcmc2 <- tagmMcmcTrain(object = d2, fcol = "xx",
                         numIter = .numIter, burnin = 0, thin = 1, numChains = 2,
                         BPPARAM = SerialParam())
  expect_equal(mcmc1, mcmc2)
  mcmc1 <- tagmMcmcProcess(mcmc1)
  mcmc2 <- tagmMcmcProcess(mcmc2)
  expect_equal(mcmc1, mcmc2)
  ans1 <- tagmPredict(object = d1, params = mcmc1, fcol = "markers")
  ans2 <- tagmPredict(object = d2, params = mcmc2, fcol = "xx")
  expect_equal(ans1, ans2, check.attributes = FALSE)
})