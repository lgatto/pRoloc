
## ## Algorithmic performance was estimated using 5-fold stratified
## ## cross-validation, which featured an additional cross-validation on
## ## each training partition in order to optimise free parameters via a
## ## grid search. This process was repeated 10 times and averaged
## ## accuracies are reported.


## ## Assess generalisation accuracy of trained classifier
## svmRegularisation_OLD<- function(object,
##                               fcol = "markers",
##                               cost = 2^(-4:4), 
##                               sigma = 10^(-2:3),
##                               times = 50,
##                               xval1 = 5,
##                               xval2 = 5,                              
##                               fun1 = mean,
##                               verbose = TRUE) {
##   ## As we are assessing generalisation, we remove
##   ## all unknown examples from the data set
##   d <- data.frame(exprs(object), markers = fData(object)[, fcol])
##   d.train <- d[d$markers != "unknown",]
##   d.train$markers <- factor(d.train$markers)
##   d.test <- d[d$markers == "unknown",]

##   .warnings <- NULL
##   .matricesL <- vector("list", length = xval1)
##   outerRes <- vector("list", length = times)

##   if (verbose) {
##     pb <- txtProgressBar(min = 0,
##                          max = xval2 * xval1 * times,
##                          style = 3)
##     ._k <- 0
##   }
##   for (.times in 1:times) {
##     .strat1 <- createFolds(d.train$markers, xval1, returnTrain = TRUE)

##     innerRes <- matrix(NA, nrow = xval1, ncol = 3)
##     colnames(innerRes) <- c("F1", "sigma", "cost")
##     for (.xval1 in 1:xval1) {    
##       .test1 <- d.train[-.strat1[[.xval1]],] ## 'unseen' test set
##       .train1 <- d.train[.strat1[[.xval1]],] ## to be used to parameter optimisation
##       .strat2 <- createFolds(.train1$markers, xval2, returnTrain = TRUE)
##       .matrixF1L <- vector("list", length=xval2)
##       for (.xval2 in 1:xval2) {
##         if (verbose) {
##           setTxtProgressBar(pb, ._k)
##           ._k <- ._k + 1
##         }
##         .train2 <- .train1[.strat2[[.xval2]], ]
##         .test2 <- .train1[-.strat2[[.xval2]], ]
##         ## The second argument in .makeF1matrix will be
##         ## used as rows, the first one for columns
##         .matrixF1 <- .makeF1matrix(list(cost = cost, sigma = sigma))
##         ## grid search for parameter selection
##         for (.cost in cost) {
##           for (.sigma in sigma) {
##             model <- svm(markers ~ ., .train2, gamma = .sigma, cost = .cost)
##             ans <- e1071:::predict.svm(model, .test2) 
##             conf <- confusionMatrix(.test2$markers, ans)$table
##             .p <- .checkNumbers(MLInterfaces:::.precision(conf))
##             .r <- .checkNumbers(MLInterfaces:::.recall(conf))
##             .f1 <- MLInterfaces:::.macroF1(.p, .r)
##             .matrixF1[as.character(.sigma), as.character(.cost)] <- .f1
##           }
##         }
##         ## we have a complete grid to be saved
##         .matrixF1L[[.xval2]] <- .matrixF1
##       }
##       ## we have xval2 grids to be summerised
##       .summaryF1 <- apply(array(do.call(cbind, .matrixF1L),
##                                 dim = c(length(sigma), length(cost), xval2)),
##                           c(1:2),
##                           fun1)
##       dimnames(.summaryF1) <- dimnames(.matrixF1)
##       .matricesL[[.xval1]] <- .summaryF1
##       .bestParams <- getBestParams(.summaryF1)[,1] ## take the first one
##       model <- svm(markers ~ ., .train1, gamma = .bestParams["sigma"], cost = .bestParams["cost"])
##       ans <- e1071:::predict.svm(model, .test1) 
##       conf <- confusionMatrix(.test1$markers, ans)$table
##       p <- .checkNumbers(MLInterfaces:::.precision(conf), tag = "precision", params = .bestParams)
##       r <- .checkNumbers(MLInterfaces:::.recall(conf), tag = "recall", params = .bestParams)
##       .f1 <- MLInterfaces:::.macroF1(p, r) ## macro F1 score for xval1 i's iteration
##       innerRes[.xval1, ] <- c(.f1, .bestParams["sigma"], .bestParams["cost"])
##     }
##     outerRes[[.times]] <- innerRes
##   }
##   if (verbose) {
##     setTxtProgressBar(pb, ._k)
##     close(pb)
##   }
  
##   .hyperparams <- list(cost = cost,
##                        sigma = sigma)
##   .design <- c(xval1 = xval1,
##                xval2 = xval2,
##                times = times)
##   ans <- new("GenRegRes",
##              hyperparameters = .hyperparams,
##              design = .design,
##              results = outerRes,
##              matrices = .matricesL)
##   if (!is.null(.warnings)) {
##     ans@log <- list(warnings = .warnings)
##     sapply(.warnings, warning)
##   }
##   return(ans)
## }

## svmPrediction_OLD <- function(object,                            
##                               assessRes,
##                               fcol = "markers") {
##   if (missing(assessRes))
##     stop("First run 'assessSvmRegularisation'")
##   ## getting best parameters
##   params <- getRegularisedParams(assessRes)
##   bkp <- NULL
##   if ("markers" != fcol) {
##     if ("markers" %in% fvarLabels(object))
##       bkp <- fData(object)$markers
##     fData(object)$markers <- fData(object)[, fcol]
##   }
##   trainInd <- which(fData(object)$markers != "unknown")  
##   ans <- MLearn(markers ~ ., t(object), svmI, trainInd,
##                 gamma = params["sigma"],
##                 cost = params["cost"])
##   fData(object)$svm <- predictions(ans)
##   fData(object)$svm.scores <- predScores(ans)
##   if (!is.null(bkp)) {
##     fData(object)$markers <- bkp
##   } else {
##     k <- which(fvarLabels(object) == "markers")
##     fData(object) <- fData(object)[, -k]
##   }
##   return(object)
## }




## plsdaRegularisation_OLD <- function(object,
##                                 fcol = "markers",
##                                 ncomp = 1:6,
##                                 times = 50,
##                                 xval1 = 5,
##                                 xval2 = 5,                              
##                                 fun1 = mean,
##                                 verbose = TRUE) {
##   ## As we are assessing generalisation, we remove
##   ## all unknown examples from the data set
##   d <- data.frame(exprs(object), markers = fData(object)[, fcol])
##   d.train <- d[d$markers != "unknown",]
##   d.train$markers <- factor(d.train$markers)
##   d.test <- d[d$markers == "unknown",]

##   .warnings <- NULL
##   .matricesL <- vector("list", length = xval1)
##   outerRes <- vector("list", length = times)
  
##   if (verbose) {
##     pb <- txtProgressBar(min = 0,
##                          max = xval2 * xval1 * times,
##                          style = 3)
##     ._k <- 0
##   }
##   for (.times in 1:times) {
##     .strat1 <- createFolds(d.train$markers, xval1, returnTrain = TRUE)
    
##     innerRes <- matrix(NA, nrow = xval1, ncol = 2)
##     colnames(innerRes) <- c("F1", "ncomp")
##     for (.xval1 in 1:xval1) {    
##       .test1 <- d.train[-.strat1[[.xval1]],] ## 'unseen' test set
##       .train1 <- d.train[.strat1[[.xval1]],] ## to be used to parameter optimisation
##       .strat2 <- createFolds(.train1$markers, xval2, returnTrain = TRUE)
##       .matrixF1L <- vector("list", length=xval2)
##       for (.xval2 in 1:xval2) {
##         if (verbose) {
##           setTxtProgressBar(pb, ._k)
##           ._k <- ._k + 1
##         }
##         .train2 <- .train1[.strat2[[.xval2]], ]
##         .test2 <- .train1[-.strat2[[.xval2]], ]
##         ## First argument used for columns
##         .matrixF1 <- .makeF1matrix(list(ncomp = ncomp))
##         ## grid search for parameter selection
##         for (.ncomp in ncomp) {
##           .x <- which(names(.train2) == fcol)
##           model <- caret::plsda(.train2[, -.x], .train2[, .x], 
##                                 probMethod = "Bayes", prior = NULL,
##                                 ncomp = .ncomp)
##           ans <- caret::predict.plsda(model, .test2[, -.x], type = "class") 
##           conf <- confusionMatrix(.test2$markers, ans)$table
##           .p <- .checkNumbers(MLInterfaces:::.precision(conf))
##           .r <- .checkNumbers(MLInterfaces:::.recall(conf))
##           .f1 <- MLInterfaces:::.macroF1(.p, .r)
##           .matrixF1[1, as.character(.ncomp)] <- .f1
##         }
##         ## we have a complete grid to be saved
##         .matrixF1L[[.xval2]] <- .matrixF1
##       }
##       ## we have xval2 grids to be summerised
##       .summaryF1 <- apply(array(do.call(cbind, .matrixF1L),
##                                 dim = c(1, length(ncomp), xval2)),
##                           c(1:2),
##                           fun1)
##       dimnames(.summaryF1) <- dimnames(.matrixF1)
##       .matricesL[[.xval1]] <- .summaryF1
##       .bestParams <- getBestParams(.summaryF1)[,1] ## take the first one      
##       .x <- which(names(.train1) == fcol)
##       model <- caret::plsda(.train1[, -.x], .train1[, .x], 
##                             probMethod = "Bayes", prior = NULL,
##                             ncomp = .ncomp)
##       ans <- caret::predict.plsda(model, .test1[, -.x], type = "class") 
##       conf <- confusionMatrix(.test1$markers, ans)$table    
##       p <- .checkNumbers(MLInterfaces:::.precision(conf), tag = "precision", params = .bestParams)
##       r <- .checkNumbers(MLInterfaces:::.recall(conf), tag = "recall", params = .bestParams)
##       .f1 <- MLInterfaces:::.macroF1(p, r) ## macro F1 score for xval1 i's iteration
##       innerRes[.xval1, ] <- c(.f1, .bestParams["ncomp"])
##     }
##     outerRes[[.times]] <- innerRes
##   }
##   if (verbose) {
##     setTxtProgressBar(pb, ._k)
##     close(pb)
##   }

##   .hyperparams <- list(ncomp = ncomp)
##   .design <- c(xval1 = xval1,
##                xval2 = xval2,
##                times = times)
##   ans <- new("GenRegRes",
##              hyperparameters = .hyperparams,
##              design = .design,
##              results = outerRes,
##              matrices = .matricesL)
##   if (!is.null(.warnings)) {
##     ans@log <- list(warnings = .warnings)
##     sapply(.warnings, warning)
##   }
##   return(ans)
## }


## plsdaPrediction_OLD <- function(object, fcol = "markers") {
##   bkp <- NULL
##   if ("markers" != fcol) {
##     if ("markers" %in% fvarLabels(object))
##       bkp <- fData(object)$markers
##     fData(object)$markers <- fData(object)[, fcol]
##   }  
##   trainInd <- which(fData(object)$markers != "unknown")
##   pls <- MLearn(markers ~ ., t(object), plsdaI, trainInd)
##   fData(object)$plsda <- predictions(pls)
##   fData(object)$plsda.scores <- predScores(pls)
##   if (!is.null(bkp)) {
##     fData(object)$markers <- bkp
##   } else {
##     k <- which(fvarLabels(object) == "markers")
##     fData(object) <- fData(object)[, -k]
##   }
##   return(object)
## }

