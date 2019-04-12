context("Annotations parameters")

## Object cc will be used in units below. Querying error-prone biomart only
## once.

library("pRolocdata")
data("hyperLOPIT2015")
hyperLOPIT2015 <- markerMSnSet(hyperLOPIT2015)
par <- setAnnotationParams(inputs = c("Mouse genes",
                                      "UniProtKB/Swiss-Prot ID"))
cc <- addGoAnnotations(hyperLOPIT2015, par,
                       namespace = "cellular_component")

test_that("Used annotation params successfully", {
  expect_true(validObject(cc))
})

context("GO annotations framework")

test_that("filtering by min or max matrix annotations gives expected results", {
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

context("clustDist function")

test_that("output from orderGoAnnotations and manually ordering clusters", {
    seed <- sample(.Machine$integer.max, 1)

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
