##' Classification parameter optimisation for the support vector
##' machine algorithm.
##'
##' Note that when performance scores precision, recall and (macro) F1
##' are calculated, any NA values are replaced by 0. This decision is
##' motivated by the fact that any class that would have either a NA
##' precision or recall would result in an NA F1 score and,
##' eventually, a NA macro F1 (i.e. mean(F1)). Replacing NAs by 0s
##' leads to F1 values of 0 and a reduced yet defined final macro F1
##' score.
##' 
##' @title svm parameter optimisation
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The feature meta-data containing marker definitions.
##' Default is \code{markers}.
##' @param cost The hyper-parameter. Default values are \code{2^-4:4}.
##' @param sigma The hyper-parameter. Default values are \code{10^(-2:3)}.
##' @param times The number of times internal cross-validation is performed.
##' Default is 100.
##' @param test.size The size of test data. Default is 0.2 (20 percent).
##' @param xval The \code{n}-cross validation. Default is 5.
##' @param fun The function used to summarise the \code{xval} macro F1 matrices.
##' @param seed The optional random number generator seed.
##' @param verbose A \code{logical} defining whether a progress bar is displayed.
##' @param ... Additional parameters passed to \code{\link{svm}} from package \code{e1071}.
##' @return An instance of class \code{"\linkS4class{GenRegRes}"}.
##' @seealso \code{\link{svmClassification}} and example therein.
##' @aliases svmRegularisation svmOptimization
##' @author Laurent Gatto
svmOptimisation <- function(object,
                            fcol = "markers",
                            cost = 2^(-4:4), 
                            sigma = 10^(-3:2),
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
  colnames(.results) <- c("F1", "sigma", "cost") 
  
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
      .matrixF1 <- makeF1matrix(list(cost = cost, sigma = sigma))
      ## grid search for parameter selection
      for (.cost in cost) {
        for (.sigma in sigma) {
          model <- svm(markers ~ ., .train2,
                       gamma = .sigma, cost = .cost, ...)
          ans <- predict(model, .test2)
          conf <- confusionMatrix(ans, .test2$markers)$table
          .p <- checkNumbers(MLInterfaces:::.precision(conf))
          .r <- checkNumbers(MLInterfaces:::.recall(conf))
          .f1 <- MLInterfaces:::.macroF1(.p, .r, naAs0 = TRUE)
          .matrixF1[as.character(.sigma), as.character(.cost)] <- .f1
        }
      }
      ## we have a complete grid to be saved
      .matrixF1L[[.xval]] <- .matrixF1
    }
    ## we have xval grids to be summerised
    .summaryF1 <- summariseMatList(.matrixF1L, fun)
    .f1Matrices[[.times]] <- .summaryF1
    .bestParams <- getBestParams(.summaryF1)[1:nparams, 1] ## takes a random best param
    model <- svm(markers ~ ., .train1,
                 gamma = .bestParams["sigma"],
                 cost = .bestParams["cost"], ...)
    ans <- predict(model, .test1) 
    .cmMatrices[[.times]] <- conf <- confusionMatrix(ans, .test1$markers)$table ## NEW
    p <- checkNumbers(MLInterfaces:::.precision(conf),
                      tag = "precision", params = .bestParams)
    r <- checkNumbers(MLInterfaces:::.recall(conf),
                      tag = "recall", params = .bestParams)
    f1 <- MLInterfaces:::.macroF1(p, r, naAs0 = TRUE) ## macro F1 score for .time's iteration
    .results[.times, ] <- c(f1, .bestParams["sigma"], .bestParams["cost"])
  }
  if (verbose) {
    setTxtProgressBar(pb, ._k)
    close(pb)
  }
  
  .hyperparams <- list(cost = cost,
                       sigma = sigma)
  .design <- c(xval = xval,
               test.size = test.size,
               times = times)

  ans <- new("GenRegRes",
             algorithm = "svm",
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

svmOptimization <-
  svmOptimisation

svmRegularisation <- function(...) {
  .Deprecated(msg = "This function has been replaced by 'svmOptimisation'.")
  svmOptimisation(...)  
}

##' Classification using the support vector machine algorithm.
##'
##' @title svm classification
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param assessRes An instance of class
##'     \code{"\linkS4class{GenRegRes}"}, as generated by
##'     \code{\link{svmOptimisation}}.
##' @param scores One of \code{"prediction"}, \code{"all"} or
##'     \code{"none"} to report the score for the predicted class
##'     only, for all classes or none.
##' @param cost If \code{assessRes} is missing, a \code{cost} must be
##'     provided.
##' @param sigma If \code{assessRes} is missing, a \code{sigma} must
##'     be provided.
##' @param fcol The feature meta-data containing marker definitions.
##'     Default is \code{markers}.
##' @param ... Additional parameters passed to \code{\link{svm}} from
##'     package \code{e1071}.
##' @return An instance of class \code{"\linkS4class{MSnSet}"} with
##'     \code{svm} and \code{svm.scores} feature variables storing the
##'     classification results and scores respectively.
##' @author Laurent Gatto
##' @aliases svmPrediction
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' ## reducing parameter search space and iterations 
##' params <- svmOptimisation(dunkley2006, cost = 2^seq(-2,2,2), sigma = 10^seq(-1, 1, 1),  times = 3)
##' params
##' plot(params)
##' f1Count(params)
##' levelPlot(params)
##' getParams(params)
##' res <- svmClassification(dunkley2006, params)
##' getPredictions(res, fcol = "svm")
##' getPredictions(res, fcol = "svm", t = 0.75)
##' plot2D(res, fcol = "svm")
svmClassification <- function(object,                            
                              assessRes,
                              scores = c("prediction", "all", "none"),
                              cost,
                              sigma,
                              fcol = "markers",
                              ...) {
  scores <- match.arg(scores)
  if (missing(assessRes)) {
    if (missing(cost) | missing(sigma))
      stop("First run 'svmOptimisation' or set 'cost' and 'sigma' manually.")
    params <- c("cost" = cost,
                "sigma" = sigma)
  } else {
    params <- getParams(assessRes)
    if (is.na(params["cost"]))
      stop("No 'cost' found.")
    if (is.na(params["sigma"]))
      stop("No 'sigma' found.")
  }
  trainInd <- which(fData(object)[, fcol] != "unknown")
  form <- as.formula(paste0(fcol, " ~ ."))
  ans <- MLearn(form, t(object), svmI, trainInd,
                gamma = params["sigma"],
                cost = params["cost"],
                ...)
  ## Would be better to check of these columns exist
  fData(object)$svm <- predictions(ans)
  if (scores == "all") {
    scoreMat <- predScores(ans)
    colnames(scoreMat) <- paste0(colnames(scoreMat), ".svm.scores")
    fData(object)$svm.all.scores <- scoreMat
  } else if (scores == "prediction") {
    fData(object)$svm.scores <- predScore(ans)
  } ## else scores is "none"
  object@processingData@processing <- c(processingData(object)@processing,
                                        paste0("Performed svm prediction (", 
                                               paste(paste(names(params), params, sep = "="),
                                                     collapse = " "), ") ",
                                               date()))
  if (validObject(object))
    return(object)
}

svmPrediction <- function(...) {
  .Deprecated(msg = "This function has been replaced by 'svmClassification'.")
  svmClassification(...)
}

