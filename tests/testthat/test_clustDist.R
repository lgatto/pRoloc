context("clustDist")

test_that("output from orderGoAnnotations and manually ordering clusters", {
  library("pRolocdata")
  data("hyperLOPIT2015")
  hyperLOPIT2015 <- markerMSnSet(hyperLOPIT2015)
  par <- setAnnotationParams(inputs = c("Mus musculus", 
                                        "UniProt/Swissprot"))
  seed <- sample(.Machine$integer.max, 1)

  cc <- addGoAnnotations(hyperLOPIT2015, par, 
                     namespace = "cellular_component")
  cc <- filterMinMarkers(cc, n = 20)
  cc <- filterMaxMarkers(cc, p = .5)
  res <- orderGoAnnotations(cc, k = 1:3, p = 1/3, 
                        verbose = FALSE, seed = seed)
  orgs1 <- colnames(fData(res)$GOAnnotations)
  
  ## Manually calculate distance and re-order by min
  dd <- clustDist(cc, fcol = "GOAnnotations", k = 1:3, 
                  verbose = FALSE, seed = seed)
  minDist <- getNormDist(dd, p = 1/3)
  o <- order(minDist)
  orgs2 <- names(dd)[o]
  
  expect_equal(orgs1, orgs2)
})
  