context("Annotations parameters")

## Object cc will be used in units below. Querying error-prone biomart only
## once.

library("pRolocdata")
data("hyperLOPIT2015")
hyperLOPIT2015 <- markerMSnSet(hyperLOPIT2015)
xx <- mrkVecToMat(hyperLOPIT2015, vfcol = "markers", mfcol = "Markers")

# par <- setAnnotationParams(inputs = c("Mouse genes",
#                                       "UniProtKB/Swiss-Prot ID"))
# cc <- addGoAnnotations(hyperLOPIT2015, par,
#                        namespace = "cellular_component")
# 
# test_that("Used annotation params successfully", {
#   expect_true(validObject(cc))
# })
# 
# context("GO annotations framework")

test_that("filtering by min or max matrix annotations gives expected results", {
    ## Check that min markers specified matches output
    m1 <- 20
    cc1 <- filterMinMarkers(xx, n = m1, fcol = "Markers")
    m2 <- min(colSums(fData(cc1)$Markers))
    expect_true(m1 <= m2)

    ## Check that max markers specified matches output
    m1 <- 50
    cc2 <- filterMaxMarkers(xx, n = m1, fcol = "Markers")
    m2 <- max(colSums(fData(cc2)$Markers))
    expect_true(m1 >= m2)
})

context("clustDist function")

# test_that("output from orderGoAnnotations and manually ordering clusters", {
#     seed <- sample(.Machine$integer.max, 1)
#     
#     xx <- filterMinMarkers(xx, n = 20, fcol = "Markers")
#     xx <- filterMaxMarkers(xx, p = .5, fcol = "Markers")
#     res <- orderGoAnnotations(xx, k = 1:3, p = 1/3, fcol = "Markers",
#                               verbose = FALSE, seed = seed)
#     orgs1 <- colnames(fData(res)$GOAnnotations)
# 
#     ## Manually calculate distance and re-order by min
#     dd <- clustDist(xx, fcol = "GOAnnotations", k = 1:3,
#                     verbose = FALSE, seed = seed)
#     minDist <- getNormDist(dd, p = 1/3)
#     o <- order(minDist)
#     orgs2 <- names(dd)[o]
# 
#     expect_equal(orgs1, orgs2)
# })
