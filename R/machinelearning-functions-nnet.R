##' Classification parameter optimisation for artificial neural network
##' algorithm.
##'
##' @title nnet parameter optimisation
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The feature meta-data containing marker definitions.
##' Default is \code{markers}.
##' @param decay The hyper-parameter. Default values are \code{c(0, 10^(-1:-5))}.
##' @param size The hyper-parameter. Default values are \code{seq(1, 10, 2)}.
##' @param times The number of times internal cross-validation is performed.
##' Default is 100.
##' @param test.size The size of test data. Default is 0.2 (20 percent).
##' @param xval The \code{n}-cross validation. Default is 5.
##' @param fun The function used to summarise the \code{xval} macro F1 matrices.
##' @param seed The optional random number generator seed.
##' @param verbose A \code{logical} defining whether a progress bar is displayed.
##' @param ... Additional parameters passed to \code{\link{nnet}} from package \code{nnet}.
##' @return An instance of class \code{"\linkS4class{GenRegRes}"}.
##' @seealso \code{\link{nnetClassification}} and example therein.
##' @aliases nnetRegularisation nnetOptimization
##' @author Laurent Gatto
nnetOptimisation <- function(object,
                             fcol = "markers",
                             decay = c(0, 10^(-1:-5)),
                             size = seq(1, 10, 2),
                             times = 100,
                             test.size = .2,
                             xval = 5,                               
                             fun = mean,
                             seed,
                             verbose = TRUE,
                             ...) {

  nparams <- 2 ## 2 or 1, depending on the algorithm
  mydata <- subsetAsDataFrame(object, fcol, train = TRUE)

  if (missing(seed)) {
    seed <- sample(.Machine$integer.max, 1)
  }
  .seed <- as.integer(seed)  
  set.seed(.seed)

  ## initialise output
  .warnings <- NULL
  .f1Matrices <- vector("list", length = times) 
  .testPartitions <- .cmMatrices <- vector("list", length = times) ## NEW
  .results <- matrix(NA, nrow = times, ncol = nparams + 1)
  colnames(.results) <- c("F1", "decay", "size") 
  
  if (verbose) {
    pb <- txtProgressBar(min = 0,
                         max = xval * times,
                         style = 3)
    ._k <- 0
  }
  for (.times in 1:times) {
    .size <- ceiling(table(mydata$markers) * test.size)
    ## size is ordered according to levels, but strata
    ## expects them to be ordered as they appear in the data
    .size <- .size[unique(mydata$markers)] 
    test.idx <- strata(mydata, "markers",
                       size = .size,
                       method = "srswor")$ID_unit
    .testPartitions[[.times]] <- test.idx ## NEW
    
    .test1   <- mydata[ test.idx, ] ## 'unseen' test set
    .train1  <- mydata[-test.idx, ] ## to be used for parameter optimisation
    
    xfolds <- createFolds(.train1$markers, xval, returnTrain = TRUE)
    ## stores the xval F1 matrices
    .matrixF1L <- vector("list", length = xval)  
    
    for (.xval in 1:xval) {    
      if (verbose) {
        setTxtProgressBar(pb, ._k)
        ._k <- ._k + 1
      }
      .train2 <- .train1[ xfolds[[.xval]], ]
      .test2  <- .train1[-xfolds[[.xval]], ]    
      
      ## The second argument in makeF1matrix will be
      ## used as rows, the first one for columns
      .matrixF1 <- makeF1matrix(list(decay = decay, size = size))
      ## grid search for parameter selection
      for (.decay in decay) {
        for (.size in size) {
          model <- nnet(markers ~ ., .train2,
                       decay = .decay, size = .size,
                        trace = FALSE, ...)          
          ans <- nnet:::predict.nnet(model, .test2, type = "class")
          ## nnet:::predict.nnet does not save the factor levels,
          ## which makes confusionMatrix to fail if both arguments
          ## have different levels.
          ans <- factor(as.character(ans), levels = levels(.test2$markers))
          conf <- confusionMatrix(ans, .test2$markers)$table
          .p <- checkNumbers(MLInterfaces:::.precision(conf))
          .r <- checkNumbers(MLInterfaces:::.recall(conf))
          .f1 <- MLInterfaces:::.macroF1(.p, .r)
          .matrixF1[as.character(.size), as.character(.decay)] <- .f1
        }
      }
      ## we have a complete grid to be saved
      .matrixF1L[[.xval]] <- .matrixF1
    }
    ## we have xval grids to be summerised
    .summaryF1 <- summariseMatList(.matrixF1L, fun)
    .f1Matrices[[.times]] <- .summaryF1
    .bestParams <- getBestParams(.summaryF1)[1:nparams, 1] ## take the first one
    model <- nnet(markers ~ ., .train1,
                  decay = .bestParams["decay"],
                  size = .bestParams["size"],
                  trace = FALSE, ...)
    ans <- nnet:::predict.nnet(model, .test1, type = "class")
    ans <- factor(as.character(ans), levels = levels(.train1$markers)) ## see comment above
    .cmMatrices[[.times]] <- conf <- confusionMatrix(ans, .test1$markers)$table ## NEW    
    p <- checkNumbers(MLInterfaces:::.precision(conf),
                      tag = "precision", params = .bestParams)
    r <- checkNumbers(MLInterfaces:::.recall(conf),
                      tag = "recall", params = .bestParams)
    f1 <- MLInterfaces:::.macroF1(p, r) ## macro F1 score for .time's iteration
    .results[.times, ] <- c(f1, .bestParams["decay"], .bestParams["size"])
  }
  if (verbose) {
    setTxtProgressBar(pb, ._k)
    close(pb)
  }
  
  .hyperparams <- list(decay = decay,
                       size = size)
  .design <- c(xval = xval,
               test.size = test.size,
               times = times)

  ans <- new("GenRegRes",
             algorithm = "nnet",
             seed = .seed,
             hyperparameters = .hyperparams,
             design = .design,
             results = .results,             
             f1Matrices = .f1Matrices,
             cmMatrices = .cmMatrices, ## NEW
             testPartitions = .testPartitions, ## NEW
             datasize = list(
               "data" = dim(mydata),
               "data.markers" = table(mydata[, "markers"]),
               "train1" = dim(.train1),
               "test1" = dim(.test1),
               "train1.markers" = table(.train1[, "markers"]),
               "train2" = dim(.train2),
               "test2" = dim(.test2),
               "train2.markers" = table(.train2[, "markers"])))

  if (!is.null(.warnings)) {
    ans@log <- list(warnings = .warnings)
    sapply(.warnings, warning)
  }
  return(ans)
}

nnetOptimization <-
  nnetOptimisation

nnetRegularisation <- function(...) {
  .Deprecated(msg = "This function has been replaced by 'nnetOptimisation'.")
  nnetOptimisation(...)  
}



##' Classification using the artificial neural network
##' algorithm.
##'
##' @title nnet classification
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param assessRes An instance of class \code{"\linkS4class{GenRegRes}"},
##' as generated by \code{\link{nnetOptimisation}}.
##' @param scores One of \code{"prediction"}, \code{"all"} or \code{"none"}
##' to report the score for the predicted class only, for all cluster
##' or none.
##' @param decay If \code{assessRes} is missing, a \code{decay} must be provided.
##' @param size If \code{assessRes} is missing, a \code{size} must be provided.
##' @param fcol The feature meta-data containing marker definitions.
##' Default is \code{markers}.
##' @return An instance of class \code{"\linkS4class{MSnSet}"} with
##' \code{nnet} and \code{nnet.scores} feature variables storing the
##' classification results and scores respectively.
##' @author Laurent Gatto
##' @aliases nnetPrediction
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' ## reducing parameter search space and iterations 
##' params <- nnetOptimisation(dunkley2006, decay = 10^(c(-1, -5)), size = c(5, 10), times = 3)
##' params
##' plot(params)
##' levelPlot(params)
##' getParams(params)
##' res <- nnetPrediction(dunkley2006, params)
##' getPredictions(res, fcol = "nnet")
##' getPredictions(res, fcol = "nnet", t = 0.75)
##' plot2D(res, fcol = "nnet")
nnetClassification <- function(object,                            
                               assessRes,
                               scores = c("prediction", "all", "none"),
                               decay,
                               size,
                               fcol = "markers") {
  scores <- match.arg(scores)  
  if (missing(assessRes)) {
    if (missing(decay) | missing(size))
      stop("First run 'nnetOptimisation' or set 'decay' and 'size' manually.")
    params <- c("decay" = decay,
                "size" = size)
  } else {
    params <- getParams(assessRes)
    if (is.na(params["decay"]))
      stop("No 'decay' found.")
    if (is.na(params["size"]))
      stop("No 'size' found.")    
  }
  trainInd <- which(fData(object)[, fcol] != "unknown")  
  form <- as.formula(paste0(fcol, " ~ ."))
  ans <- MLearn(form, t(object), nnetI, trainInd,
                decay = params["decay"],
                size = params["size"])  
  fData(object)$nnet <- predictions(ans)
  if (scores == "all") {
    scoreMat <- predScores(ans)
    colnames(scoreMat) <- paste0(colnames(scoreMat), ".nnet.scores")
    fData(object) <- cbind(fData(object), scoreMat)
  } else if (scores == "prediction") {
    fData(object)$nnet.scores <- predScore(ans)
  } ## else scores is "none" 
  object@processingData@processing <-
    c(processingData(object)@processing,
      paste0("Performed nnet prediction (", 
             paste(paste(names(params), params, sep = "="),
                   collapse = " "), ") ",
             date()))
  if (validObject(object))
    return(object)
}

nnetPrediction <- function(...) {
  .Deprecated(msg = "This function has been replaced by 'nnetClassification'.")
  nnetClassification(...)
}

