context("distances")

m <- matrix(c(rep(0, 12),
              c(1, 2, 9, 10)),
            ncol = 4)
colnames(m) <- 1:4
rownames(m) <- LETTERS[1:4]

q <- matrix(c(rep(0, 4),
              rep(10, 4)),
            nrow = 2,
            byrow = TRUE)    
colnames(q) <- 1:4
rownames(q) <- LETTERS[10:11]

e <- new("MSnSet", exprs = m,
         featureData = new("AnnotatedDataFrame",
             data = data.frame(X = 1:4,
                 row.names = LETTERS[1:4]))) 

test_that("nndist errors", {
    expect_error(nndist(m, q, k = 2, dist = "euclidean"))
    expect_error(nndist(m, q, k = 2, dist = "mahalanobis"))
    expect_error(nndist(m, k = 2, dist = "foo"))
    expect_error(nndist(e, k = 2, dist = "foo"))
})

test_that("nndist matrix/msnset", {
    expect_true(validObject(e))    
    resm0 <- pRoloc:::nndist_matrix(m, k = 2, dist = "euclidean")
    resm <- nndist(m, k = 2, dist = "euclidean")
    expect_equal(resm0, resm)
    rese0 <- pRoloc:::nndist_msnset(e, k = 2, dist = "euclidean")
    rese <- nndist(e, k = 2, dist = "euclidean")
    expect_equal(rese0, rese)
    expect_equal(data.frame(resm), fData(rese)[, -1]) 
})

test_that("nndist query", {
    res <- nndist(m, q, k = 2)
    res0 <- rbind(c(1, as.numeric(dist(rbind(q[1, ], m[1, ]))),
                    2, as.numeric(dist(rbind(q[1, ], m[2, ])))),
                  c(4, as.numeric(dist(rbind(q[2, ], m[4, ]))),
                    3, as.numeric(dist(rbind(q[2, ], m[3, ])))))
    dimnames(res0) <- dimnames(res)
    expect_equal(res, res0)
})
