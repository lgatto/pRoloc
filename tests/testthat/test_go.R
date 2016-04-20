test_that("GO term/id printing", {
    expect_identical(goIdToTerm("GO:0000001"),
                     c('GO:0000001' = "mitochondrion inheritance"))
    expect_identical(goIdToTerm("GO:0000001", names = FALSE),
                     "mitochondrion inheritance")
    expect_identical(goIdToTerm(c("GO:0000001", "GO:0000002", "notvalid")),
                     c('GO:0000001' = "mitochondrion inheritance",
                       'GO:0000002' = "mitochondrial genome maintenance",
                       notvalid = NA))
    expect_identical(goTermToId("mitochondrion inheritance"),
                     c('mitochondrion inheritance' = "GO:0000001"))
    expect_identical(goTermToId("mitochondrion inheritance", names = FALSE),
                     "GO:0000001")
    expect_identical(goTermToId(c("mitochondrion inheritance", "notvalid")),
                     c('mitochondrion inheritance' = "GO:0000001", notvalid = NA))
    expect_identical(prettyGoTermId("mitochondrion inheritance"),
                     "mitochondrion inheritance (GO:0000001)")
    expect_identical(prettyGoTermId("GO:0000001"),
                     "mitochondrion inheritance (GO:0000001)")
    expect_identical(flipGoTermId("mitochondrion inheritance"),
                     c('mitochondrion inheritance' = "GO:0000001"))
    expect_identical(flipGoTermId("GO:0000001", names = FALSE),
                     "mitochondrion inheritance")
})
