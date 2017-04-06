test_that("AnnotationParams", {
    expect_error(setAnnotationParams(1),
                 "'params' must be of class 'AnnotationParams'.")
    expect_error(setAnnotationParams(inputs = 1),
                 "'inputs' must contain a species and a feature type.")
    expect_error(setAnnotationParams(inputs = LETTERS[1:2]),
                 "Couldn't find a unique species match for 'A'.")
    expect_error(setAnnotationParams(inputs = c("Homo sapiens", "foo")),
                 "Couldn't find a unique feature type match for 'foo'.")

    ap <- setAnnotationParams(inputs = c("Homo sapiens",
                                         "UniProtKB/Swiss-Prot ID"))
    expect_null(show(ap))    
    
    expect_is(ap, "AnnotationParams")
    expect_identical(ap@martname, "ENSEMBL_MART_ENSEMBL")
    expect_identical(ap@dataset, "hsapiens_gene_ensembl")
    expect_identical(ap@filter, "uniprotswissprot")
    expect_error(setAnnotationParams()) ## only in non-interactive mode

    ap2 <- getAnnotationParams()
    expect_identical(ap, ap2)
    setAnnotationParams(ap)
    ap2 <- getAnnotationParams()
    expect_identical(ap, ap2)    
})
