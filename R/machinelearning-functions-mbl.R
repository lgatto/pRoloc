mblPrediction <- function(object,
                          fcol = "markers",
                          verbose = TRUE,
                          scores = c("prediction", "all", "none"),
                          seed,
                          ...) {
  if (missing(seed)) {
    seed <- sample(.Machine$integer.max, 1)
  }
  .seed <- as.integer(seed)  
  set.seed(.seed)

  if (verbose)
    message("Using seed ", .seed, ".")
  
  scores <- match.arg(scores)
  
  levs <- unique(as.character(fData(object)[, fcol]))
  levs <- levs[levs != "unknown"]
  regs <- preds <- vector("list", length = length(levs))
  names(regs) <- names(preds) <- levs
  
  ans <- lapply(levs, function(l) {
    if (verbose)
      message("Analysing binary ", l, " cluster...")
    obj <- object
    mrk <- as.character(fData(obj)[, fcol])
    mrk[mrk != "unknown" & mrk != l] <- paste0("not-", l)
    fData(obj)[, fcol] <- factor(mrk)
    ## if we want the same partitions for each binary classification,
    ## seed = .seed should be passed to svmRegularisation
    regs[[l]] <<- reg <- svmRegularisation(obj, fcol,
                                           seed = .seed,
                                           verbose = verbose, ...)
    preds[[l]] <<- pred <- svmPrediction(obj, reg, fcol = fcol,
                                         scores = "all")
    res <- fData(pred)[, grep("scores", fvarLabels(pred))]
    res
  })  
  ans <- do.call(cbind, ans)
  fData(object) <- cbind(fData(object), ans)
  stopifnot(validObject(object))
  
  ## preparing new MSnSet with mlb scores
  .eset <- as.matrix(fData(object)[, grep("scores", fvarLabels(object))])
  .fd <- featureData(object)[, -grep("scores", fvarLabels(object))]
  object2 <- new("MSnSet",
                 exprs = .eset,
                 featureData = .fd,
                 processingData = processingData(object))
  stopifnot(validObject(object2))
  
  ## running svm with new scores
  if (verbose)
    message("Prediction on combined binary results...")
  reg2 <- svmRegularisation(object2, fcol = fcol, verbose = verbose, ...)
  object2 <- svmPrediction(object2, reg2, fcol = fcol, scores = scores)
  object2@processingData@processing <- c(processingData(object2)@processing,
                                         paste0("Performed multiple binary learnig (mbl) prediction ",
                                                date()))
  if (validObject(object2))
    return(object2)
}
