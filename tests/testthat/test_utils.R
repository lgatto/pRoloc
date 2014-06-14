context("utils")

test_that("subsetAsDataFrame preserves col/rownames", {
    library("MSnbase")
    set.seed(1)
    data(itraqdata)
    itraqdata <- itraqdata[1:10, ]
    msnset <- quantify(itraqdata ,method = "trap",
                       reporters = iTRAQ4, verbose = FALSE)
    fData(msnset)$markers <- sample(LETTERS[1:5],
                                    nrow(itraqdata),
                                    replace=TRUE)
    .fcol <- "markers"
    dfr <- pRoloc:::subsetAsDataFrame(msnset, fcol = .fcol)
    expect_equal(featureNames(msnset), rownames(dfr))
    expect_equal(c(sampleNames(msnset), .fcol), colnames(dfr))
})


test_that("getPredictions with different thresholds", {
    library("pRolocdata")
    data(tan2009r1)
    ## using parameters as estimated in vignette
    ## (although these are from markers, not pd.markers)
    res <- svmClassification(dunkley2006, fcol = "pd.markers", sigma = 0.1, cost = 0.5)
    classres <- as.character(fData(res)$svm)
    allpreds <- getPredictions(res, fcol = "svm", t = 0, verbose = FALSE)
    expect_equal(classres, allpreds)
    preds <- getPredictions(res, fcol = "svm", t = 1, verbose = FALSE)
    mrkrs <- getMarkers(res, "pd.markers", verbose = FALSE)
    names(mrkrs) <- NULL
    expect_equal(preds, mrkrs)
    tserr <- ts <- tapply(fData(res)$svm.scores, fData(res)$svm, median)
    names(tserr)[1] <- "wrongnames"
    expect_error(getPredictions(res, fcol = "svm", t = tserr, verbose = FALSE))
    preds <- getPredictions(res, fcol = "svm", t = ts, verbose = FALSE)
    for (k in unique(classres)) {
        ## all thresholded predictions must be > that class score
        i <- preds == k
        expect_true(all(fData(res)[i, "svm.scores"] >= ts[k]))
        ## predicted but not in thresholded predictions must have
        ## lower scores than class scores
        j <- allpreds == k
        expect_true(all(fData(res)[j & !i, "svm.scores"] < ts[k]))
    }        
})



## Class FeatComp, function compfnames and helper functions .compStrings and 
## .calcCompNumbers

test_that("FeatComp", {
    ## tests for class FeatComp
    expect_error(pRoloc:::.FeatComp(name = 1, common = "a", unique1 = "b", 
                           unique2 = "c", all = TRUE))
    expect_error(pRoloc:::.FeatComp(name = "test", common = 1, unique1 = "b", 
                           unique2 = "c", all = TRUE))
    expect_error(pRoloc:::.FeatComp(name = "test", common = "a", unique1 = 1, 
                           unique2 = "c", all = TRUE))
    expect_error(pRoloc:::.FeatComp(name = "test", common = "a", unique1 = "b", 
                           unique2 = 1, all = TRUE))
    expect_error(pRoloc:::.FeatComp(name = "test", common = "a", unique1 = "b", 
                           unique2 = "c", all = "d"))
    test <- pRoloc:::.FeatComp(name = "test", common = "a", 
                    unique1 = letters[2:4], unique2 = "c", all = TRUE)
    expect_is(test, "FeatComp")
    expect_equal(test@name, "test")
    expect_equal(test@common, "a")
    expect_equal(test@unique1, c("b", "c", "d"))
    expect_equal(test@unique2, "c")
    expect_true(test@all)
})

test_that("FeatComp_.compStrings", {
    ## create vectors with simple strings
    string1 <- c("a", "b", "c", "d")
    string2 <- c("c", "d", "e", "f")
    string3 <- c("f", "e")
    expect_equal(length(pRoloc:::.compStrings(string1, string2)), 1)
    ## test slots
    expect_true(pRoloc:::.compStrings(string1, string2)@all, TRUE)
    expect_false(pRoloc:::.compStrings(string1, string2, all = FALSE)@all)
    expect_equal(pRoloc:::.compStrings(string1, string2)@name, "all")
    expect_equal(pRoloc:::.compStrings(string1, string2, 
                                all = FALSE, name = "ut")@name, "ut")
    ## test common features
    expect_equal(pRoloc:::.compStrings(string1, string2)@common, c("c", "d"))
    expect_equal(pRoloc:::.compStrings(string1, string3)@common, character())
    expect_equal(pRoloc:::.compStrings(string2, string3)@common, c("e", "f"))
    ## test unique features for argument string1
    expect_equal(pRoloc:::.compStrings(string1, string2)@unique1, c("a", "b"))
    expect_equal(pRoloc:::.compStrings(string2, string3)@unique1, c("c", "d"))
    expect_equal(pRoloc:::.compStrings(string1, string3)@unique1, 
                                c("a", "b", "c", "d"))
    ## test unique features for argument string2
    expect_equal(pRoloc:::.compStrings(string1, string2)@unique2, c("e", "f"))
    expect_equal(pRoloc:::.compStrings(string2, string3)@unique2, character())
    expect_equal(pRoloc:::.compStrings(string1, string3)@unique2, c("f", "e"))
})

test_that("FeatComp_.calcCompNumbers", {
    ## test .calcCompNumbers which creates matrix 
    ## create objects for input for function .calcCompNumbers
    flist11 <- pRoloc:::.compStrings(string1, string1, FALSE, "ut11")
    flist12 <- pRoloc:::.compStrings(string1, string2, FALSE, "ut12")
    flist13 <- pRoloc:::.compStrings(string1, string3, FALSE, "ut13")
    flist23 <- pRoloc:::.compStrings(string2, string3, FALSE, "ut23")
    flist <- list(flist11, flist12, flist13, flist23)
    ## start testing
    expect_true(is.matrix(pRoloc:::.calcCompNumbers(flist)))
    expect_equal(dim(pRoloc:::.calcCompNumbers(flist)), c(4, 3))
    expect_equal(rownames(pRoloc:::.calcCompNumbers(flist)), 
                            c("ut11", "ut12", "ut13", "ut23"))
    expect_equal(colnames(pRoloc:::.calcCompNumbers(flist)), 
                            c("common", "unique1", "unique2"))
    expect_equal(as.vector(pRoloc:::.calcCompNumbers(flist)[1, ]), c(4, 0, 0)) 
    expect_equal(as.vector(pRoloc:::.calcCompNumbers(flist)[2, ]), c(2, 2, 2))
    expect_equal(as.vector(pRoloc:::.calcCompNumbers(flist)[3, ]), c(0, 4, 2))
    expect_equal(as.vector(pRoloc:::.calcCompNumbers(flist)[4, ]), c(2, 2, 0))
})

test_that("FeatComp_compfnames", {
    ## test function compfnames
    ## create three test MSnSets
    obj1 <- new("MSnSet",
                exprs = matrix(1, ncol = 4, nrow = 4, 
                            dimnames = list(letters[1:4])),             
                featureData = new("AnnotatedDataFrame",
                                data = data.frame(markers = LETTERS[1:4],
                                        row.names = letters[1:4])))
    obj2 <- new("MSnSet", 
                exprs = matrix(1, ncol = 4, nrow = 4, 
                            dimnames = list(letters[3:6])),
                featureData = new("AnnotatedDataFrame",
                                data = data.frame(markers = LETTERS[3:6],
                                        row.names = letters[3:6])))
    obj3 <- new("MSnSet", 
                exprs = matrix(1, ncol = 4, nrow = 4, 
                            dimnames = list(letters[5:8])),
                featureData = new("AnnotatedDataFrame",
                                data = data.frame(markers = rep(LETTERS[6:7],2),
                                        row.names = letters[5:8])))
    ## start testing function compfnames
    expect_error(compfnames(flist, obj2, NULL, NULL))
    expect_error(compfnames(obj1, flist, NULL, NULL))
    expect_error(compfnames(obj1, obj2, "pd.markers", "markers"))
    expect_error(compfnames(obj1, obj2, "markers", "pd.markers"))
    expect_error(compfnames(obj1, obj2, "markers", "markers", "verbose"))
    ## create object which compares featureNames of obj1 and obj2
    comp12 <- compfnames(obj1, obj2, NULL, "markers", FALSE)
    expect_equal(length(comp12), 1)
    expect_equal(comp12[[1]]@name, "all")
    expect_equal(comp12[[1]]@common, c("c", "d"))
    expect_equal(comp12[[1]]@unique1, c("a", "b"))
    expect_equal(comp12[[1]]@unique2, c("e", "f"))
    expect_true(comp12[[1]]@all)
    comp12 <- compfnames(obj1, obj2, "markers", verbose = FALSE)
    expect_equal(comp12[[4]]@common, "c")
    expect_equal(comp12[[4]]@unique1, character())
    expect_equal(comp12[[4]]@unique2, character())
    expect_equal(comp12[[5]]@common, "d")
    ## create object which compares featureNames of obj1 and obj3
    comp13 <- compfnames(obj1, obj3, "markers", "markers", FALSE)
    expect_equal(length(comp13), 7)
    expect_equal(comp13[[1]]@common, character())
    expect_equal(comp13[[1]]@unique1, c("a", "b", "c", "d"))
    expect_equal(comp13[[1]]@unique2, c("e", "f", "g", "h"))
    expect_equal(comp13[[2]]@common, character())
    expect_equal(comp13[[2]]@unique1, "a")
    expect_equal(comp13[[2]]@unique2, character())
    expect_equal(comp13[[7]]@common, character())
    expect_equal(comp13[[7]]@unique1, character())
    expect_equal(comp13[[7]]@unique2, c("f", "h"))
    ## create object whcih compares featureNames of obj2 and obj3
    comp23 <- compfnames(obj2, obj3, "markers", "markers", FALSE)
    expect_equal(length(comp23), 6)
    expect_equal(comp23[[1]]@common, c("e", "f"))
    expect_equal(comp23[[1]]@unique1, c("c", "d"))
    expect_equal(comp23[[1]]@unique2, c("g", "h"))
    expect_equal(comp23[[2]]@common, character())
    expect_equal(comp23[[2]]@unique1, "c")
    expect_equal(comp23[[2]]@unique2, character())
    expect_equal(comp23[[6]]@common, character())
    expect_equal(comp23[[6]]@unique1, character())
    expect_equal(comp23[[6]]@unique2, c("f", "h"))                        
})

