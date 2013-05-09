context("Data consistence")

data(dunkley2006)

test_that("subsetting", {
  mrk <- markerSet(dunkley2006)
  nm <- sum(fData(dunkley2006)$markers == "unknown")
  expect_equal(nrow(dunkley2006) - nm, nrow(mrk))
})
