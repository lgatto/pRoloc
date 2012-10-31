##' Hyper-parameters regularisation for the partial least square
##' distcriminant analysis algorithm.
##'
##' @title plsda regularisation
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The feature meta-data containing marker definitions.
##' Default is \code{markers}.
##' @param ncomp The hyper-parameter. Default values are \code{1:6}.
##' @param times The number of times internal cross-validation is performed.
##' Default is 100.
##' @param test.size The size of test data. Default is 0.2 (20 percent).
##' @param xval The \code{n}-cross validation. Default is 5.
##' @param fun The function used to summarise the \code{times} macro F1 matrices.
##' @param seed The optional random number generator seed.
##' @param verbose A \code{logical} defining whether a progress bar is displayed.
##' @param ... Additional parameters passed to \code{\link{plsda}} from package \code{caret}.
##' @return An instance of class \code{"\linkS4class{GenRegRes}"}.
##' @seealso \code{\link{plsdaPrediction}} and example therein.
##' @author Laurent Gatto
plsdaRegularisation <- function(object,
                                fcol = "markers",
                                ncomp = 1:6,
                                times = 100,
                                test.size = .2,
                                xval = 5,                               
                                fun = mean,
                                seed,
                                verbose = TRUE,
                                ...) {

  nparams <- 1 ## 2 or 1, depending on the algorithm
  mydata <- subsetAsDataFrame(object, fcol, train = TRUE)

  if (missing(seed)) {
    seed <- sample(.Machine$integer.max, 1)
  }
  .seed <- as.integer(seed)  
  set.seed(.seed)

  ## initialise output
  .warnings <- NULL
  .matrices <- vector("list", length = times) 
  .results <- matrix(NA, nrow = times, ncol = nparams + 1)
  colnames(.results) <- c("F1", "ncomp") 
  
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
        .matrixF1 <- makeF1matrix(list(ncomp = ncomp))
        ## grid search for parameter selection
          for (.ncomp in ncomp) {
          .x <- which(names(.train2) == "markers")
          model <- caret::plsda(.train2[, -.x], .train2[, .x], 
                                probMethod = "Bayes", prior = NULL,
                                ncomp = .ncomp, ...)
          ans <- caret::predict.plsda(model, .test2[, -.x], type = "class") 
          conf <- confusionMatrix(ans, .test2$markers)$table
          .p <- checkNumbers(MLInterfaces:::.precision(conf))
          .r <- checkNumbers(MLInterfaces:::.recall(conf))
          .f1 <- MLInterfaces:::.macroF1(.p, .r)
          .matrixF1[1, as.character(.ncomp)] <- .f1
        }
        ## we have a complete grid to be saved
        .matrixF1L[[.xval]] <- .matrixF1
      }
    ## we have xval grids to be summerised
    .summaryF1 <- summariseMatList(.matrixF1L, fun)
    .matrices[[.times]] <- .summaryF1
    .bestParams <- getBestParams(.summaryF1)[1:nparams, 1] ## take the first one
    .x <- which(names(.train1) == "markers")
    model <- caret::plsda(.train1[, -.x], .train1[, .x], 
                          probMethod = "Bayes", prior = NULL,
                          ncomp = .bestParams["ncomp"], ...)
    ans <- caret::predict.plsda(model, .test1[, -.x], type = "class") 
    conf <- confusionMatrix(ans, .test1$markers)$table
    p <- checkNumbers(MLInterfaces:::.precision(conf),
                      tag = "precision", params = .bestParams)
    r <- checkNumbers(MLInterfaces:::.recall(conf),
                      tag = "recall", params = .bestParams)
    f1 <- MLInterfaces:::.macroF1(p, r) ## macro F1 score for .time's iteration
    .results[.times, ] <- c(f1, .bestParams["ncomp"])
  }
  if (verbose) {
    setTxtProgressBar(pb, ._k)
    close(pb)
  }
  
  .hyperparams <- list(ncomp = ncomp)
  .design <- c(xval = xval,
               test.size = test.size,
               times = times)

  ans <- new("GenRegRes",
             algorithm = "plsda",
             seed = .seed,
             hyperparameters = .hyperparams,
             design = .design,
             results = .results,
             matrices = .matrices,
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


##' Classification using the partial least square
##' distcriminant analysis algorithm.
##'
##' @title plsda prediction
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param assessRes An instance of class \code{"\linkS4class{GenRegRes}"},
##' as generated by \code{\link{plsdaRegularisation}}.
##' @param scores One of \code{"prediction"}, \code{"all"} or \code{"none"}
##' to report the score for the predicted class only, for all cluster
##' or none.
##' @param ncomp If \code{assessRes} is missing, a \code{ncomp} must be provided.
##' @param fcol The feature meta-data containing marker definitions.
##' Default is \code{markers}.
##' @return An instance of class \code{"\linkS4class{MSnSet}"} with
##' \code{plsda} and \code{plsda.scores} feature variables storing the
##' classification results and scores respectively.
##' @author Laurent Gatto
##' @examples
##' \dontrun{
##' ## not running this one for time considerations
##' library(pRolocdata)
##' data(dunkley2006)
##' ## reducing parameter search space and interations 
##' reg <- plsdaRegularisation(dunkley2006, ncomp = c(3, 10),  times = 2)
##' reg
##' plot(reg)
##' levelPlot(reg)
##' getRegularisedParams(reg)
##' res <- plsdaPrediction(dunkley2006, reg)
##' getPredictions(res, fcol = "plsda")
##' getPredictions(res, fcol = "plsda", t = 0.75)
##' plot2D(res, fcol = "plsda")
##' }
plsdaPrediction <- function(object,
                            assessRes,
                            scores = c("prediction", "all", "none"),
                            ncomp,
                            fcol = "markers") {
  scores <- match.arg(scores)  
  if (missing(assessRes)) {
    if (missing(ncomp))
      stop("First run 'plsdaRegularisation' or set 'ncomp' manually.")
    params <- c("ncomp" = ncomp)
  } else {
    params <- getRegularisedParams(assessRes)
    if (is.na(params["ncomp"]))
      stop("No 'ncomp' found.")

  }
  trainInd <- which(fData(object)[, fcol] != "unknown")
  form <- as.formula(paste0(fcol, " ~ ."))
  ans <- MLearn(form, t(object), plsdaI, trainInd,
                ncomp = params["ncomp"])
  fData(object)$plsda <- predictions(ans)
  if (scores == "all") {
    scoreMat <- predScores(ans)
    colnames(scoreMat) <- paste0(colnames(scoreMat), ".plsda.scores")
    fData(object) <- cbind(fData(object), scoreMat)
  } else if (scores == "prediction") {
    fData(object)$plsda.scores <- predScore(ans)
  } ## else scores is "none" 
  object@processingData@processing <- c(processingData(object)@processing,
                                        paste0("Performed pls-da prediction (", 
                                               paste(paste(names(params), params, sep = "="),
                                                     collapse = " "), ") ",
                                               date()))
  if (validObject(object))
    return(object)
}

