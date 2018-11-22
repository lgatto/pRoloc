context("Tagm MAP tests")

data(dunkley2006)
d2 <- d1 <- dunkley2006
i <- which(fvarLabels(d1) == "markers")
stopifnot(length(i) == 1)
fvarLabels(d2)[i] <- "xx"

.times <- 2
.seed <- 1


test_that("TAGM consistency", {
  .numIter <- 2
  map1 <- tagmMapTrain(object = d1, fcol = "markers",
                       numIter = .numIter, seed = .seed)
  map2 <- tagmMapTrain(object = d2, fcol = "xx",
                       numIter = .numIter, seed = .seed)
  expect_equal(map1, map2)
  ans1 <- tagmPredict(object = d1, params = map1, fcol = "markers")
  ans2 <- tagmPredict(object = d2, params = map2, fcol = "xx")
  expect_equal(ans1, ans2, check.attributes = FALSE)
})