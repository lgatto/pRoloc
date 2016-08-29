context("plot2D using visualTest")
## reference figures are in inst/figs/.
library("visualTest")

test_that("Negative control", {
    expect_false(
        isSimilar(system.file("figs/plot2D-null-black.png", package = "pRoloc"),
                  system.file("figs/plot2D-null.png", package = "pRoloc")))
    expect_false(
        isSimilar(system.file("figs/plot2D-null-black.png", package = "pRoloc"),
                  system.file("figs/plot2D.png", package = "pRoloc")))
})

test_that("Positive control", {
    expect_true(
        isSimilar(system.file("figs/plot2D-null.png", package = "pRoloc"),
                  system.file("figs/plot2D-null.png", package = "pRoloc")))
    expect_true(
        isSimilar(system.file("figs/plot2D.png", package = "pRoloc"),
                  system.file("figs/plot2D.png", package = "pRoloc")))
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL)", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL)
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-null.png"))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, col = 'black')", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, col = "black")
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-null-black.png"))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015)", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015)
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D.png"))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, cex = 3)", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, cex = 3)
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-null-cex.png"))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, pch = 3)", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, pch = 3)
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-null-pch.png"))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, pch = 19, col = 'black')", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, pch = 19, col = "black")
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-null-pch-black.png"))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, method = 'scree')", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, method = "scree")
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-scree.png"))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, method = 'hexbin')", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, method = "hexbin")
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-hexbin.png"))
    unlink(tmpf)
})

test_that("plot2D(hyperLOPIT2015, method = 'lda')", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, method = "lda")
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-lda.png"))
    unlink(tmpf)
})
