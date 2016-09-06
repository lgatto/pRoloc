context("plot2D using visualTest")
## reference figures are in inst/figs/.
library("visualTest")

## see https://github.com/MangoTheCat/visualTest/issues/19
thr <- ifelse(Sys.info()[["nodename"]] == "elyacin", 1e-3, 25)

test_that("Negative control", {
    expect_false(
        isSimilar(system.file("figs/plot2D-null-black.png", package = "pRoloc"),
                  system.file("figs/plot2D-null.png", package = "pRoloc"),
                  threshold = thr))
    expect_false(
        isSimilar(system.file("figs/plot2D-null-black.png", package = "pRoloc"),
                  system.file("figs/plot2D.png", package = "pRoloc"),
                  threshold = thr))
})

test_that("Positive control", {
    expect_true(
        isSimilar(system.file("figs/plot2D-null.png", package = "pRoloc"),
                  system.file("figs/plot2D-null.png", package = "pRoloc"),
                  threshold = thr))
    expect_true(
        isSimilar(system.file("figs/plot2D.png", package = "pRoloc"),
                  system.file("figs/plot2D.png", package = "pRoloc"),
                  threshold = thr))
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL)", {
    ref <- system.file("figs/plot2D-null.png", package = "pRoloc")
    data(hyperLOPIT2015)   
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL)
    dev.off()
    expect_true(isSimilar(tmpf, ref, threshold = thr))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, col = 'black')", {
    ref <- system.file("figs/plot2D-null-black.png", package = "pRoloc")
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, col = "black")
    dev.off()
    expect_true(isSimilar(tmpf, ref, threshold = thr))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015)", {
    ref <- system.file("figs/plot2D.png", package = "pRoloc")
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015)
    dev.off()
    expect_true(isSimilar(tmpf, ref))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, cex = 3)", {
    ref <- system.file("figs/plot2D-null-cex.png", package = "pRoloc")
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, cex = 3)
    dev.off()
    expect_true(isSimilar(tmpf, ref, threshold = thr))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, pch = 3)", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, pch = 3)
    dev.off()
    expect_true(isSimilar(tmpf,
                          system.file("figs/plot2D-null-pch.png", package = "pRoloc"),
                          threshold = thr))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, pch = 19, col = 'black')", {
    ref <- system.file("figs/plot2D-null-pch-black.png", package = "pRoloc")
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, pch = 19, col = "black")
    dev.off()
    expect_true(isSimilar(tmpf, ref, threshold = thr))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, method = 'scree')", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, method = "scree")
    dev.off()
    expect_true(isSimilar(tmpf,
                          system.file("figs/plot2D-scree.png", package = "pRoloc"),
                          threshold = thr))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, method = 'hexbin')", {
    ref <- system.file("figs/plot2D-hexbin.png", package = "pRoloc")
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, method = "hexbin")
    dev.off()
    expect_true(isSimilar(tmpf, ref, threshold = thr))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, method = 'lda')", {
    ref <- system.file("figs/plot2D-lda.png", package = "pRoloc")
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, method = "lda")
    dev.off()
    expect_true(isSimilar(tmpf, ref, threshold = thr))
    unlink(tmpf)
})
