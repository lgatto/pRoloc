context("GO marker annotation tests")

test_that("output from orderGoMarkers and the ouput from manually ordering GO markers using clustDist are the same", {
  library("pRolocdata")
  data("hyperLOPIT2015")
  hyperLOPIT2015 <- markerMSnSet(hyperLOPIT2015)
  
  cc <- addGoMarkers(hyperLOPIT2015, params, 
                     namespace = "cellular_component")
  cc <- filterMinMarkers(cc)
  cc <- filterMaxMarkers(cc)
  
  seed <- sample(.Machine$integer.max, 1)
  
  ## method 1 (using orderGoMarkers)
  res <- orderGoMarkers(cc, k = 1:3, p = 1/3, 
                        verbose = FALSE, seed = seed)
  resGO1 <- colnames(fData(res)$GOMarkers)
  
  ## method 2 (manually calculating and re-ordering)
  dd <- clustDist(cc, fcol = "GOMarkers", k = 1:3, 
                  verbose = FALSE, seed = seed)
  minDist <- getNormDist(dd, p = 1/3)
  o <- order(minDist)
  fData(cc)$GOMarkers <- fData(cc)$GOMarkers[, o]
  resGO2 <- colnames(fData(cc)$GOMarkers)
  
  expect_equal(resGO1, names(minDist))
  expect_equal(resGO1, resGO2)
  
})