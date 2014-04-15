context("organelle markers")


test_that("adding markers", {
    set.seed(1)
    x <- new("MSnSet",
             exprs = matrix(1, ncol = 4, nrow = 10,
                 dimnames = list(1:10)),             
             featureData = new("AnnotatedDataFrame",
                 data = data.frame(
                     markers = rep(c("A", "B"), each = 5),
                     id = c(LETTERS[1:7], "A", "Z", "Z"),
                     id2 = c(letters[1:7], "A", "Z", "Z"))))    
    m <- sample(paste0("ORG", 1:5), 26, replace = TRUE)
    names(m) <- LETTERS    
    ## Error: Detected an existing 'markers' feature column.
    expect_error(addMarkers(x, n)) 
    fData(x)$markers <- NULL
    ## Error: No markers found. Are you sure that the feature names match?
    expect_error(addMarkers(x, n))     
    ## Case 1 (default): using feature names (all match)
    featureNames(x2) <- LETTERS[1:10]
    fData(x2)$markers <- NULL
    x2 <- addMarkers(x2, m)
    expect_true(all(m[1:10] == fData(x2)$markers))
    ## Case 2: using an fcol 
    x3 <- addMarkers(x, m, fcol = "id")
    expect_true(all(m[as.character(fData(x3)$id)] ==
                    fData(x3)$markers))
    ## Case 3: using fcol, with unknowns
    x4 <- addMarkers(x, m, fcol = "id2")
    expect_true(all(c(rep("unknown", 7), as.character(fData(x3)$markers[8:10]))
                    == fData(x4)$markers))
    ## Case 4: markers not present in featureNames/fcol
    names(m)[1:10] <- letters[1:10]
    x5 <- addMarkers(x, m, fcol = "id")
    expect_true(all(c(rep("unknown", 8), as.character(fData(x3)$markers[9:10]))
                    == fData(x5)$markers))        
})
