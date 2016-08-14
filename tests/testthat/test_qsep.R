context("QSep infrastructure")

test_that("QSep disance matrices", {
    library("pRolocdata")
    data(dunkley2006)
    fData(dunkley2006)$markers0 <- fData(dunkley2006)$markers
    fData(dunkley2006)$markers <- "unknown"
    match("ER", fData(dunkley2006)$markers0)
    er <- grep("ER", fData(dunkley2006)$markers0)
    fData(dunkley2006)[er, "markers"] <- 
        fData(dunkley2006)[er, "markers0"]
    qs <- QSep(dunkley2006)
    erl <- grep("ER lumen", fData(dunkley2006)$markers0)
    erm <- grep("ER membrane", fData(dunkley2006)$markers0)
    expect_identical(qsep(qs, FALSE)[1, 1],
                     mean(dist(exprs(dunkley2006)[erl, ])))
    expect_identical(qsep(qs, FALSE)[2, 2],
                     mean(dist(exprs(dunkley2006)[erm, ])))
    tmp <- as.matrix(dist(exprs(dunkley2006)[er, ]))
    diag(tmp) <- NA
    expect_identical(qsep(qs, FALSE)[2, 1],
                     mean(tmp[featureNames(dunkley2006)[erm],
                              featureNames(dunkley2006)[erl]]))
    expect_identical(qsep(qs, FALSE)[1, 2],
                     mean(tmp[featureNames(dunkley2006)[erm],
                              featureNames(dunkley2006)[erl]]))
    expect_identical(qsep(qs)[1, 2],
                     mean(tmp[featureNames(dunkley2006)[erm],
                              featureNames(dunkley2006)[erl]]) /
                     mean(dist(exprs(dunkley2006)[erl, ])))
    expect_identical(qsep(qs)[2, 1],
                     mean(tmp[featureNames(dunkley2006)[erm],
                              featureNames(dunkley2006)[erl]]) /
                     mean(dist(exprs(dunkley2006)[erm, ])))
})
