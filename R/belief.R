diffscores <- function(object, scols = "svm.scores", diffcol) {
  stopifnot(inherits(object, "MSnSet"))
  sidx <- grep(scols, fvarLabels(object))
  if (length(sidx) == 0)
    stop("No score 'svm.scores' columns found.")
  if (length(sidx) == 1) {
    msg <- c("\nOnly 1 ", fvarLabels(object)[sidx] , " column found.\n",
             "Are you sure you specified 'scores = \"all\"' when ",
             "running you classification?")
    stop(msg)
  }
  smat <- fData(object)[, sidx]
  ans <- apply(smat, 1, function(x) {
    .prob <- sort(x, decreasing = TRUE)[1:2]
    (.prob[1] - .prob[2])/.prob[1]    
  })
  if (missing(diffcol))
    diffcol <- sub("scores", "diff", scols)
  fData(x)[, diffcol] <- ans
  if (validObject(x))
    return(x)
}
