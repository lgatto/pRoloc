

## Separator for datasets
## 
## @title separateDataSet
## @param theDataSet: the table (item x variables)
## @param fcol : if exists, the name of the column in theDataSet corresponding to the labels
## @return resultSet : a dataset (matrix) and the labels (vector)
separateDataSet <- function(theDataSet, fcol = NULL){
  .theData <- NULL
  .theLabels <- NULL
  
  if (is.null(fcol) == FALSE)
    {
    ind <-which(colnames(theDataSet) != fcol)
    .theData <- theDataSet[,ind]
    .theLabels <- theDataSet[,-ind]
  }
  else{
    .theData <- theDataSet
  }
  
  return(list(theData = .theData, theLabels = .theLabels))
}


## Constructor for datasets
## 
## @title constructDataSet
## @param theDataSet: the table (item x variables)
## @param theLabel (optional) : the label corresponding to items
## @return resultSet : a dataset
constructDataSet <- function(theDataSet, theLabel = NULL){
  nCol <- dim(theDataSet)[2]
  nRow <- dim(theDataSet)[1]
  theDataSet <- matrix(unlist(theDataSet), ncol=nCol, byrow=F)
  
  resultSet <- NULL;
  resultSet$var <- theDataSet;
  resultSet$nbInd <- nRow;
  if (is.null(theLabel)==F){
    resultSet$label <- theLabel 
    resultSet$liLabel <- unique(resultSet$label);
    resultSet$nbLabel <- length(resultSet$liLabel); 
  }
  return(resultSet);	
}
 


## Control the validity of the matrix inversion method 
## combined with the regularisation method.
## 
## @title controlParameters
## @param inv The inversion method
## @param reg The type of regularisation
## @return A list of (reg, inv). Those are unchanged if the combinaison is implemented
## or set to default values (inv="Inversion Cholesky", reg = "tikhonov")
## @author Samuel Wieczorek
controlParameters <- function(inv, reg) {
  ## LG: inv and reg are arg.match'ed in the calling function.
  ## Verification de la bonne strategie d'inversion/regularisation.
  strategies <- c("Inversion Cholesky/tikhonov",
                  "Inversion Cholesky/none",
                  "Moore Penrose/tikhonov",
                  "Moore Penrose/none",
                  "solve/tikhonov",
                  "solve/none",
                  "svd/tikhonov",
                  "svd/none",
                  "svd/trunc")
  comb <- paste(inv, reg, sep = "/")
  if (!pmatch(comb, strategies, nomatch = FALSE)) {
    msg <- paste0(comb, " method of matrix inversion/regularisation is not implemented.\n",                  
                  "Using default 'Inversion Cholesky/tikhonov' regularisation\n")
    warning(msg)
    inv <- "Inversion Cholesky"
    reg <- "tikhonov"
  }  
  return(list(inv=inv, reg=reg))
}

## Learn a model from a dataset
## 
## @title learnOneClass
## @param oneClassLearningSet: the training examples to process
## @param pSigma: The variance of the kernel used to approximate the Laplace-Beltrami operator
## @param pRegul: a parameter to tune the importance of the regularization during matrix inversion
## @return A matrix which is the inverse of the Laplace-Beltrami proxy
## @author Thomas Burger, xxxxxx
learnOneClass <- function(oneClassLearningSet, pSigma, inv, reg, pRegul) {
  learningSetSize <- dim(oneClassLearningSet)[1]
  rbfkernel <- rbfdot(1/pSigma^2)  
  if (inv == "Inversion Cholesky" && reg == "tikhonov") { # INV : Inverse de Cholesky
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    kell <- kell + pRegul * diag(1, learningSetSize, learningSetSize) # REG Tikhonov regularisation
    cholkell <- chol(kell) ## pd here with the regularisation\n
    ans <- chol2inv(cholkell)
  } else if (inv == "Inversion Cholesky" && reg == "none") { # INV : svd
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    cholkell <- chol(kell) ## pd here with the regularisation\n
    ans <- chol2inv(cholkell)
  } else if (inv == "Moore Penrose" && reg == "tikhonov") { # INV : Moore-Penrose    
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    kell <- kell + pRegul * diag(1, learningSetSize, learningSetSize) # REG Tikhonov regularisation
    ans <- ginv(kell) ## ginv calculates the Moore-Penrose generalized inverse of a matrix\n
  } else if (inv == "Moore Penrose" && reg == "none"){ # INV : svd
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    ans <- ginv(kell) ## ginv calculates the Moore-Penrose generalized inverse of a matrix\n
  } else if (inv == "solve" && reg == "tikhonov") {
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    kell <- kell + pRegul * diag(1, learningSetSize, learningSetSize) # REG Tikhonov regularisation
    ans <- solve(kell)
  } else if (inv == "solve" && reg == "none") {
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    ans <- solve(kell)
  } else if (inv == "svd" && reg == "tikhonov") { # INV : svd
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    kell <- kell + pRegul * diag(1, learningSetSize, learningSetSize) # REG Tikhonov regularisation
    s <- svd(kell)
    D <- diag(1/s$d)
    U <- s$u
    ans <- U %*% D %*% t(U) ## was Kinv
  } else if (inv == "svd" && reg == "none") { # INV : svd
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    s <- svd(kell)
    D <- diag(1/s$d)
    U <- s$u
    ans <- U %*% D %*% t(U) ## was Kinv
  } else if (inv == "svd" && reg == "trunc") { # INV : svd
    ## cat("the pLambda parameter is in ]0;1]. It represents the % of the total eigenvalues kept.")\n
    kell <- kernelMatrix(rbfkernel, oneClassLearningSet)
    s <- svd(kell)
    ss <- cumsum(s$d)
    len <- length(ss)
    ss <- ss/ss[len]
    rk <- min(which(pmin((ss - pRegul), (rep(0, len))) == 0))
    D <- as.matrix(diag(1/s$d)[1:rk, 1:rk])
    U <- s$u[, 1:rk]
    ans <- U %*% D %*% t(U) ## was Kinv
  }
  return(ans)
}


## Gives the value k(x)^t K^{-1} k(x) for a set of points x
## 
## @title testStep
## @param testSet: a dataset structure
## @param trainingResults: a structure such as provided by trainingStep(), that contains K^{-1}
## @param loop: a string that indicates the type of implementation for the loop
## Available values are :'LOOP_R', 'MATRIX_R', 'LOOP_CXX' and 'BLOCKS_R'. Default is 'LOOP_CXX'
## @param nbBlocks: number of blocks of data to decompose the calculus in the BLOCKS_R implementation
## @return setOfAllPerturbations: the "distance" of each point x to each class
## @author Thomas Burger, Laetitia Chapel, Samuel Wieczorek
testStep <- function(testSet, trainingResults, loop = "LOOP_CXX", nbBlocks = 3) {
  #browser()
  setOfAllPerturbations <- matrix(rep(0, trainingResults$nbLabel * testSet$nbInd),
                                  ncol = trainingResults$nbLabel)
	if (loop =="LOOP_R")
  {
    for (iClass in 1:trainingResults$nbLabel) {
      ## We begin with the class that have the smaller sigma value
      i <- sort(trainingResults$listOfSigma, index.return = TRUE)$ix[iClass]
      rbfkernel <- rbfdot(1/trainingResults$listOfSigma[i]^2)
      kTestPoints <- kernelMatrix(rbfkernel, trainingResults$trExPerCl[[i]], testSet$var);  	 
      for (k in 1:testSet$nbInd){
        setOfAllPerturbations[k, iClass] <-
          t(kTestPoints[,k]) %*% trainingResults$listOfLearnedModels[[i]] %*% (kTestPoints[,k])
      }
    } 
    
  }
  else if  (loop =="MATRIX_R"){
    for (iClass in 1:trainingResults$nbLabel) {
      ## We begin with the class that have the smaller sigma value
      i <- sort(trainingResults$listOfSigma, index.return = TRUE)$ix[iClass]
      rbfkernel <- rbfdot(1/trainingResults$listOfSigma[i ]^2)
      kTestPoints <- kernelMatrix(rbfkernel, trainingResults$trExPerCl[[i]], testSet$var);     
      ## In the perTurbo paper [ref], the perturbation is calculated as 1-p
      ## and each sample is associated to the class with the least perturbation
      ## Here, we record only the 'p' value in order to have one operation less
      ## and each sample will be associated to the class with the bigger perturbation value
       setOfAllPerturbations[, iClass] <-
          diag(t(kTestPoints) %*% trainingResults$listOfLearnedModels[[i]] %*% kTestPoints)
      }
  }
  else if (loop =="LOOP_CXX"){
    for (iClass in 1:trainingResults$nbLabel) {
      ## We begin with the class that have the smaller sigma value
      i <- sort(trainingResults$listOfSigma, index.return = TRUE)$ix[iClass]
      rbfkernel <- rbfdot(1/trainingResults$listOfSigma[i]^2)
      
      kTestPoints <- kernelMatrix(rbfkernel, trainingResults$trExPerCl[[i]], testSet$var);  	 
      ## In the perTurbo paper [ref], the perturbation is calculated as 1-p
      ## and each sample is associated to the class with the least perturbation
      ## Here, we record only the 'p' value in order to have one operation less
      ## and each sample will be associated to the class with the bigger perturbation value
      setOfAllPerturbations[, iClass] <- loopInTestStep(kTestPoints,trainingResults$listOfLearnedModels[[i]],testSet$nbInd)
          }
  }
  else if (loop =="BLOCKS_R"){
    if (nbBlocks > testSet$nbInd){
      print ("The number of decomposition blocks must be greater or equal to the number of individuals in the dataset.")
      stop()
    }
    
    for (iClass in 1:trainingResults$nbLabel) {
      ## We begin with the class that have the smaller sigma value
      i <- sort(trainingResults$listOfSigma, index.return = TRUE)$ix[iClass]
      rbfkernel <- rbfdot(1/trainingResults$listOfSigma[i]^2)
      
      #nbBlocks <- 3
      sizeBlock <- ceiling(testSet$nbInd /nbBlocks)
    
      for (k in 1:nbBlocks){
        if (k == nbBlocks){
          rangeblock <- seq(from = (k-1)*sizeBlock + 1, to = testSet$nbInd)
          }else{
        rangeblock <- seq(from = (k-1)*sizeBlock + 1, to = k*sizeBlock)
          }
        #utile pour reformer une matrice avec les dimensions compatibles avec trainingResults$trExPerCl[[i]]
        #sinon, testSet$var[rangeblock,] renvoie un vecteur-colonne
        if (length(rangeblock) == 1)
          {
          kTestPoints <- kernelMatrix(rbfkernel, trainingResults$trExPerCl[[i]], as.matrix(t(testSet$var[rangeblock,])))
          }else{
          kTestPoints <- kernelMatrix(rbfkernel, trainingResults$trExPerCl[[i]], testSet$var[rangeblock,])
          }
         
        setOfAllPerturbations[rangeblock,iClass] <- diag(t(kTestPoints) %*% trainingResults$listOfLearnedModels[[i]] %*% kTestPoints);
        }
    }
  }
 
  return(setOfAllPerturbations)
}


estimateVariance <- function(learningSet, nbInd) {
  ## Idle function by now.
  ## Could be used in the future to adapt the parameters of a model to the training examples
  ## 
  ## rule of thumb number 1 
  kN <- round(log(nbInd)+1);
  ## problem here: matrix may be too big: can not allocate... error 
  distMat <- as.matrix(dist(learningSet, method="euclidean", diag=TRUE, upper=TRUE));
  distMat[which(distMat == 0)] <- max(distMat);
  sortedDist <- apply(distMat,1,sort);
  keptDist <- sortedDist[1:kN,];
  expSigma <-  mean(keptDist);
  ## ------------------------
  return(expSigma)
}


## Training step in which all the models are learned interatively 
## using function learnOneClass(). The variance of the kernel used to approximate the Laplace-Beltrami operator is 
## computed via estimateVariance() function
## 
## @title trainingStepAutoSigma
## @param learningSet: a dataset
## @param pRegul A parameter to tune the importance of the regularization during matrix inversion
## @return trainingResults : a structure containing several parts:
##      - trainingResults@nbLabel: nb of labels per class
##    	- trainingResults@listOfLearnedModels: K^{-1} computed for each class
##    	- trainingResults@trExPerCl: list of points per class
##    	- trainingResults@listOfSigma: list of Sigma used for each class
## @author Thomas Burger, xxxxxx
trainingStepAutoSigma <- function(theLearningSet, inv, reg, pRegul) {
  trainingResults <- NULL
  trainingResults$listOfLearnedModels <- list()
  trainingResults$trExPerCl <- list()
  trainingResults$listOfSigma <- rep(0, theLearningSet$nbLabel) #vector's initialisation\n
  trainingResults$nbLabel <- theLearningSet$nbLabel #usefull for testStep();\n

  for (i in 1:theLearningSet$nbLabel) {
    thisClass <- theLearningSet$var[which(theLearningSet$label == theLearningSet$liLabel[i]), ]
    trainingResults$listOfSigma[i] <- estimateVariance(thisClass, dim(thisClass)[1])
    trainingResults$trExPerCl <- c(trainingResults$trExPerCl, list(thisClass))
  }

  ## In the case where we want only one sigma value\n
  vSigma <- min(trainingResults$listOfSigma)
  for (i in 1:theLearningSet$nbLabel) {
    ## only one sigma value\n
    trainingResults$listOfSigma[i] <- vSigma
    ## thisClassModel <- learnOneClass(trainingResults@trExPerCl[[i]], vSigma , pRegul, inv)
    thisClassModel <- learnOneClass(trainingResults$trExPerCl[[i]],
                                    trainingResults$listOfSigma[i], 
                                    inv, reg, pRegul)    
    trainingResults$listOfLearnedModels <-
      c(trainingResults$listOfLearnedModels, list(thisClassModel))
  }
  return(trainingResults)
}

## Training step in which all the models are learned interatively
## using function learnOneClass(). Id but here we use sigma as an input instead of
## being computed inside the loop
## 
## @title trainingStep
## @param learningSet: a dataset
## @param pRegul A parameter to tune the importance of the regularization during matrix inversion
## @param pSigma The variance of the kernel used to approximate the Laplace-Beltrami operator
## @return trainingResults : a structure containing several parts:
##      - trainingResults@nbLabel: nb of labels per class
##      - trainingResults@listOfLearnedModels: K^{-1} computed for each class
##    	- trainingResults@trExPerCl: list of points per class
##    	- trainingResults@listOfSigma: list of Sigma used for each class
## @author Thomas Burger, xxxxxx
trainingStep <- function(theLearningSet, pSigma, inv, reg, pRegul) {
  trainingResults <- NULL
  trainingResults$listOfLearnedModels <- list()
  trainingResults$trExPerCl <- list()
  trainingResults$listOfSigma <- rep(pSigma, theLearningSet$nbLabel)
  trainingResults$nbLabel <- theLearningSet$nbLabel ## usefull for testStep()

  for (i in 1:theLearningSet$nbLabel) {
    thisClass <- theLearningSet$var[which(theLearningSet$label == theLearningSet$liLabel[i]), ]
    thisClassModel <- learnOneClass(thisClass, trainingResults$listOfSigma[i], inv, reg, pRegul)
    trainingResults$trExPerCl <- c(trainingResults$trExPerCl, list(thisClass))
    trainingResults$listOfLearnedModels <- c(trainingResults$listOfLearnedModels, list(thisClassModel))
  }
  return(trainingResults)
}


## Given a set of parameters, train perTurbo and 
## assess the accuracy rate on a training set.
## 
## @title trainingPerTurbo
## @param theDataSet: the table (num[1: nb points, 1: nb var]) (mandatory)
## @param theLabel : the label corresponding to the examples (mandatory)
## @param inv: method used to inverse the K matrix (mandatory)
## @param reg : type of regularisation method (mandatory)
## @param pRegul parameter for the regularisation (mandatory)
## @param propTraining: proportion of the training set to use for training,
## the other part is for testing. DEFAULT: 0.7. (optional)
## @param sigma:"AUTO" if we pick the "best" sigma, vector of sigma otherwise. Default: AUTO (optional)
## @param scaled: if the data are in [0,1]. DEFAULT: TRUE (optional)
## @return trainingResults : a structure containing several parts:
##      - trainingResults@nbLabel: nb of labels per class
##      - trainingResults@listOfLearnedModels: K^{-1} computed for each class
##      - trainingResults@trExPerCl: list of points per class
##    	- trainingResults@listOfSigma: list of Sigma used for each class
## @author Thomas Burger, xxxxxx
trainingPerTurbo <- function(markers, train2, sigma , inv, reg, pRegul, scaled=FALSE) {
  
  learningSet <- constructDataSet(train2, markers)
 if (scaled == TRUE) {
    learningSet$var <- scale.default(learningSet$var,
                       center = apply(learningSet$var, 2, min),
                       scale = apply(learningSet$var, 2, max, na.rm = TRUE) - apply(learningSet$var, 2, min, na.rm = TRUE))
    }
  
  ## Training with fixed sigma or a sigma given by the user
  if (sigma == "AUTO") {
    trModel <- trainingStepAutoSigma(learningSet,  inv, reg, pRegul)
  } else {
    trModel <- trainingStep(learningSet, sigma, inv, reg, pRegul)
  }
  return(trModel)
}


## Given a set of parameters, test the model given by perTurbo on a test set.
## 
## @title testPerTurbo
## @param testSet: the table (num[1: nb points, 1: nb var]) (mandatory)
## @param trModel : the model learnt by TrainingPerTurbo
## @param markers: labels of testSet
## @return list of estimated classes for the test data
## @author Thomas Burger, xxxxxx
testPerTurbo <- function(trModel, markers, testSet) {
  testSet <- constructDataSet(testSet, markers)
  TestResult <- testStep(testSet, trModel, loop="LOOP_R")
  ListOfEstimatedClasses <- testSet$liLabel[apply(TestResult, 1, which.max)]
  return(ListOfEstimatedClasses)
}


## Given a set of parameters, test the model given by perTurbo on a test set.
## 
## @title predictionPerTurbo
## @param trModel : the model learnt by TrainingPerTurbo
## @param testSet: the table (num[1: nb points, 1: nb var]) (mandatory)
## @param markers: labels of testSet
## @return xxxxx
## @author Thomas Burger, xxxxxx
predictionPerTurbo <- function(trModel, markers, preTestSet) {
  testSet <- constructDataSet(preTestSet, markers)
  TestResult <- testStep(testSet, trModel)
  return(TestResult)
}
