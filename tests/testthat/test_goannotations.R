context("go annotations framework")

test_that("filtering by min or max matrix annotations gives expected results", {
  library("pRolocdata")
  data("hyperLOPIT2015")
  hyperLOPIT2015 <- markerMSnSet(hyperLOPIT2015)
  par <- setAnnotationParams(inputs = c("Mus musculus", 
                                        "UniProtKB/Swiss-Prot ID"))
  cc <- addGoAnnotations(hyperLOPIT2015, par, 
                     namespace = "cellular_component")
  
  ## Check that min markers specified matches output
  m1 <- 20
  cc1 <- filterMinMarkers(cc, n = m1)
  m2 <- min(colSums(fData(cc1)$GOAnnotations))
  expect_true(m1 <= m2)
  
  ## Check that max markers specified matches output
  m1 <- 50
  cc2 <- filterMaxMarkers(cc, n = m1)
  m2 <- max(colSums(fData(cc2)$GOAnnotations))
  expect_true(m1 >= m2)
  
})
