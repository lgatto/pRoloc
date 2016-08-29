context("plot2D using visualTest")
## reference figures are in inst/figs/.

test_that("Negative control", {
    expect_false(
        isSimilar("~/dev/00_github/pRoloc/inst/figs/plot2D-null-black.png",
                  "~/dev/00_github/pRoloc/inst/figs/plot2D-null.png"))
    expect_false(
        isSimilar("~/dev/00_github/pRoloc/inst/figs/plot2D-null-black.png",
                  "~/dev/00_github/pRoloc/inst/figs/plot2D.png"))
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL)", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL)
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-null.png"))
})

test_that("plot2D(hyperLOPIT2015, fcol = NULL, col = 'black')", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015, fcol = NULL, col = "black")
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D-null-black.png"))
})

test_that("plot2D(hyperLOPIT2015)", {
    data(hyperLOPIT2015)
    tmpf <- paste0(tempfile(), ".png")
    png(tmpf)
    plot2D(hyperLOPIT2015)
    dev.off()
    expect_true(isSimilar(tmpf, "~/dev/00_github/pRoloc/inst/figs/plot2D.png"))
})

