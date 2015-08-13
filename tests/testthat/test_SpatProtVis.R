context("SpatProtVis class and constructor")

test_that("Consistency with plot2D", {
    library("pRolocdata")
    data(dunkley2006)
    ## in the interest of time, don't use t-SNE
    m <- c("PCA", "MDS", "kpca")
    vis <- SpatProtVis(dunkley2006, methods = m)
    show(vis)
    plot(vis)
    dev.off()
    x <- plot2D(dunkley2006, method = "kpca", plot = FALSE)
    expect_equal(x, vis@vismats[["kpca"]])
    x <- plot2D(dunkley2006, method = "PCA", plot = FALSE)
    expect_equal(x, vis@vismats[["PCA"]])
    x <- plot2D(dunkley2006, method = "MDS", plot = FALSE)
    expect_equal(x, vis@vismats[["MDS"]])

    ## Setting method arguments
    margs <- c(list(kpar = list(sigma = 0.1)),
               list(kpar = list(sigma = 1.0)),
               list(kpar = list(sigma = 10)),
               list(kpar = list(sigma = 100)))
    vis <- SpatProtVis(dunkley2006,
                       methods = rep("kpca", 4),
                       methargs = margs)
    x <- plot2D(dunkley2006, method = "kpca", plot = FALSE,
                methargs = list(kpar = list(sigma = 0.1)))
    expect_equal(x, vis@vismats[[1]])
    x <- plot2D(dunkley2006, method = "kpca", plot = FALSE,
                methargs = list(kpar = list(sigma = 1)))
    expect_equal(x, vis@vismats[[2]])
    x <- plot2D(dunkley2006, method = "kpca", plot = FALSE,
                methargs = list(kpar = list(sigma = 10)))
    expect_equal(x, vis@vismats[[3]])
    x <- plot2D(dunkley2006, method = "kpca", plot = FALSE,
                methargs = list(kpar = list(sigma = 100)))
    expect_equal(x, vis@vismats[[4]])
                
    ## Multiple PCA plots but different PCs
    dims <- list(c(1, 2), c(3, 4))
    vis <- SpatProtVis(dunkley2006, methods = c("PCA", "PCA"), dims = dims)
    x <- plot2D(dunkley2006, method = "PCA", dims = 1:2, plot = FALSE)
    expect_equal(x, vis@vismats[[1]])
    x <- plot2D(dunkley2006, method = "PCA", dims = 3:4, plot = FALSE)
    expect_equal(x, vis@vismats[[2]])
})
