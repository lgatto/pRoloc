## Function to get new phenotypes
## INPUTS: object = output from phenoDiscoTracking, groupSize = default is 5 
## which is used to define minimum number of proteins per new phenotype
## jc = correlation cooefficient, currently no other value but 1 is supported
getNewClusters <- function(object, 
                           groupSize = 5,
                           jc = 1) {  
  X <- object[,1]$originalX
  namesX <- table(unlist(apply(object, 2, function(z) rownames(z$X))))
  namesX <- names(which(namesX == dim(object)[2]))
  index <- unlist(object[1,])
  history <- sapply(namesX,
                          function(z) as.vector(index[which(names(index) == z)]), 
                          simplify = TRUE)
  history <- t(history)
  rownames(history) <- namesX
  corMatrix <- simil(history, history, method="eJaccard")
  ## Here we look at the correlation between profiles
  ## We start look at the 1st protein & identify which proteins are highly 
  ## correlated with in i.e. that has a correlation coefficent > jc
  ## We then group these proteins together to form a new phenotype
  if (jc == 1) {
    getNames <- apply(corMatrix, 1, function(z) names(z[z==1]))	
  } else {
    getNames <- apply(corMatrix, 1, function(z) names(z[z > jc])) 
  } 
  if (length(getNames) > 0) {
    if (class(getNames) == "list") {
      getNamesUnique <- unique(getNames)
      group <- getNamesUnique[unlist(lapply(getNamesUnique, 
                                            function(z) length(z) >= groupSize))]
      coords <- lapply(group, function(x) 
                             t(sapply(x, function(z) X[rownames(X) == z,])))
    } else {
      group <- NULL
      coords <- NULL
    }
  } else {
    coords <- NULL
  }
  list(coords = coords, 
       protIDs = group)
}

## Function for performing outlier detection - L. Breckels - 16/06/2011  
## L: labelled, X: unlabelled, N: number of iterations, p: significance level 
gmmOutlier <- function(L, X, N = 500, p=0.05) {
  ## Make X Matrix
  X <- matrix(X, ncol=ncol(L), dimnames=list( c(rownames(X)), c(colnames(X))))
  ## Generate Null
  ## Need justification of options for selection G here
  ## Re-test
  if (nrow(L) < 30) {
    if (nrow(L) < 10) {
      gmm0<- Mclust(L, G=1) 
    } else {
      gmm0<- Mclust(L, G=1:3)
    } 
  } else {
    gmm0<- Mclust(L)
  }	
  if(gmm0$G==1) { 
    mat <- mahalanobis(X, gmm0$parameters$mean,gmm0$parameters$variance$sigma[,,1]) 		
    ## If the cluster number in the data is 1 use the Mahalanobis distance
  } else {
    W <- WN <- a<- vector()
    for (i in 1:N) {
      s<-which(rmultinom(1, size=1, prob=(gmm0$parameters$pro))==1)
      NP<-rmultnorm(1, mu=gmm0$parameters$mean[,s], 
                    vmat=gmm0$parameters$variance$sigma[,,s]) 		
      ## Generate new profile (NP) from the data
      es<-do.call("estep", c(list(data=rbind(NP, L)), gmm0))   					
      ## ELSE use the estep of the EM algorithm to determine model parameters
      W[i]<-(-2*(es$loglik - gmm0$loglik)) 
      ## Generate the test statistic, W, for round N (build up a distribution 
      ## of W over N rounds)
    }
  }
  ## Test unlabelled
  ## Test for G>1
  if (gmm0$G!=1) {	
    WA<-W[order(W)][round((1-p)*length(W))] ## Determine W alpha 	
    for (i in (1:nrow(X)) ) {
      esN<-do.call("estep", c(list(data=rbind(L,X[i, ])), gmm0)) 	
      ## Use the estep of the EM algorithm to determine model 
      ## parameters for the unlabelled profile
      WN[i]<- (-2*(esN$loglik - gmm0$loglik)) 
      ## Calculate the test statistic, W, for the unlabelled profile
    }
    TF <- WN > WA
    TF <- replace(TF, TF=="TRUE", "outlier")
    TF <- replace(TF, TF=="FALSE", "classified")	
    fullList <- data.frame(X, TF)
    outliers <- fullList[which(fullList[,3]=="outlier"), 1:2]
    outliers <- as.matrix(outliers)
    classified <- fullList[which(fullList[,3]=="classified"), 1:2]
    classified <- as.matrix(classified)    
    if (length(outliers)==0) {outliers <- NULL}
    if (length(classified)==0) {classified <- NULL}
    list(outliers=outliers, classified=classified,fullList=fullList,gmm0=gmm0) 
    ## Compare WN with WA to determine if outlier/class member
  } else {	## Test for G=1
    chi <- qchisq(df=ncol(L)-1, 1-p)
    TF <- mat > chi
    TF <- replace(TF, TF=="TRUE", "outlier")
    TF <- replace(TF, TF=="FALSE", "classified")
    fullList <- data.frame(X, TF)
    outliers <- fullList[which(fullList[,3]=="outlier"), 1:2]
    outliers <- as.matrix(outliers)
    classified <- fullList[which(fullList[,3]=="classified"), 1:2]
    classified <- as.matrix(classified)    
    if (length(outliers) == 0) {
      outliers <- NULL
    }
    if (length(classified) == 0) {
      classified <- NULL
    }	
    list(outliers = outliers, 
         classified = classified, 
         fullList = fullList, 
         gmm0 = gmm0)
  }
}

## tracking.R (L Breckels - Date last modified 15/08/2012)
## 
## This is a core part of the phenoDisco algorihtm. This function (1) transforms 
## the data by PCA, then loops over k classes and (2) clusters each set k U X 
## using a GMM (NB: we should probably change to hierarchical so as not to confuse 
## the GMMs used in the outlier detection part) at this stage of clustering we the
## cluster number each profile assigned to is noted (2B) for later use in defining
## phenotypes and also candidates for outlier detection are identified (2A), 
## following candidate identification the 'gmmOutlier' function is called and 
## candidates are merged or rejected accordingly.
##
tracking <- function(data, alpha = 0.05, markerCol = "markers") {
  indexCol <- which(colnames(fData(data))==markerCol)
  ## ===STAGE 1=== TRANSFORM DATA by PCA to reduce data complexity
  k <- names(table(fData(data)[,indexCol]))
  k <- k[which(k!="unknown")]
  k <- sample(k) ## Sample k to avoid bias (order of classes *does* affect cluster 
                 ## members but should not affect the ID of new phenotypes)
  pca <- prcomp(exprs(data), scale=T)$x[,1:2]
  tmp <- lapply(k, function(z) fData(data)[,indexCol]==z)
  L <- lapply(tmp, function(z) pca[z , ])
  X <- pca[fData(data)[,indexCol]=="unknown", ]
  originalL <- L
  originalX <- X
  names(originalL) <- names(L) <- k

  ## Define data structure
  track <- newClassData <- outlierList <- classifiedList <- gmmCluster <- list()
  percentL <-  vector()
  ## Run over k known organelle k
  for (i in 1:length(k)) {
    outlierDetect <- list()
    ## ===STAGE 2=== Cluster the set (L(k) U X) using a GMM and -
    ## (A)  Identify potential new members/candidates of organelle class k
    ## (B)  Get cluster IDs for defining phenotypes at a later stage  
    kUx <- rbind(L[[i]], X)
    gmmCluster <- Mclust(kUx, modelNames = c("EEE","EEV","VEV","VVV"))
    track[[i]] <- gmmCluster$classification # Get cluster number for each prot
    Nclass <- nrow(L[[i]])
    classifyL <- track[[i]][1:Nclass]  # Cluster number for each L
    classifyX <- track[[i]][(Nclass+1):length(track[[i]])] # and for each X
    ## Do not include clusters where there is only 1 member! Can not call
    ## outlier detection when cluster contains only one profile! Set min group
    ## size to 5 members. 
    clusterID <- names(which(table(classifyL) >= 5))
    percentL <- sapply(clusterID, function(z) 
                       (sum(classifyL == z)/sum(track[[i]] == z)*100))
    ## If one of the clusters contains ONLY labelled members amd no unassigned  
    ## proteins we remove the cluster as outlier detection can not be called!
    if (sum(percentL==100) > 0) {clusterID <- names(which(percentL!=100))}
    ## An if/else statement is here for the instance when we have no
    ## profiles to consider i.e. we get 100% labelled profiles in *ALL* 
    ## clusters in L(k) U X
    if (length(clusterID) < 1) {
      newClassData[[i]] <- L[[i]]; X <- X;
      outlierList[[i]] <- 0; classifiedList[[i]] <- 0
    } else {
      ## Get candidates for outlier detection i.e. unlabelled profiles that are
      ## found in the same clusters as the labelled proteins.  
      indexCandidates <- unlist(lapply(clusterID, function(z) 
                                       which(classifyX==z)))
      candidates <- X[indexCandidates,]
      if (!is.matrix(candidates)) { # Might be able to remove this.....
        candidates <- matrix(candidates, ncol=ncol(L[[i]]), 
                             dimnames = list(c(
                               rownames(X)[indexCandidates]),
                               c(colnames(X))))
      }
      updateX <- X[-indexCandidates, ]
      ## ===Stage 3=== Call outlier detection on candidates
      outlierDetect <- gmmOutlier(L[[i]], candidates, p=alpha)
      ## ===Stage 4=== Merge/reject candidates. Profiles classified merged
      ## with the current class and rejected members are returned to X
      X <- rbind(updateX, outlierDetect$outliers)
      newClassData[[i]] <- rbind(L[[i]], outlierDetect$classified)
    }
  }
  for (i in 1:length(L)) {
    L[[i]] <- newClassData[[i]]
  }
  
  list(trackHistory = track, 
       X = X,
       originalX = originalX,
       L = L,
       originalL = originalL,
       k = k)
}

## Convenience function for updating MSnSet with new phenotypes
updateobject  <- function(MSnSetToUpdate,
                          newPhenotypes,
                          newClasses,
                          originalMarkerColumnName = "markers",
                          oldMarkerColumnName = "markers",
                          newMarkerColumnName = "newMarkers") {
  newobject <- MSnSetToUpdate
  indexOld <- which(colnames(fData(newobject)) == oldMarkerColumnName)
  indexNew <- ncol(fData(newobject)) + 1
  indexOriginal <- which(colnames(fData(newobject)) == originalMarkerColumnName)
  
  ## Now need to add new newPhenotypes and delete these proteins from unlabelled
  if (length(newPhenotypes$protIDs) > 0) {
    if (class(newPhenotypes$coords)!="NULL") {
      indOrigM <- length(table(fData(newobject)[,indexOriginal]))
      indOldM <- length(table(fData(newobject)[,indexOld]))
      
      nInd <- indOldM - indOrigM
      newLabels <- sapply(1:length(newPhenotypes$coords), 
                          function(z) paste("Phenotype", nInd+z)) 
      names(newPhenotypes$coords) <- newLabels
      names(newPhenotypes$protIDs) <- newLabels
      
      index <- lapply(newPhenotypes$protIDs, function(z) 
                      as.vector(sapply(z, function(x) which(featureNames(newobject)==x))))
      
      
      fData(newobject)[,indexNew] <- 
        as.character(fData(newobject)[,indexOld])
      
      for (i in 1:length(newPhenotypes$protIDs)) {
        fData(newobject)[index[[i]], indexNew] <- 
          rep(x=names(index)[i], times=length(index[[i]]))
      }
      
      index <- lapply(newClasses, function(z) 
                      as.vector(sapply(z, function(x) which(featureNames(newobject)==x))))
      
      
      for (i in 1:length(newClasses)) {
        fData(newobject)[index[[i]], indexNew] <-
          rep(x=names(index)[i], times=length(index[[i]]))
      }
      
      fData(newobject)[,indexNew] <- 
        as.factor(fData(newobject)[,indexNew])
      colnames(fData(newobject))[indexNew] <- newMarkerColumnName
    }
  } else {
    
    fData(newobject)[,indexNew] <- 
      as.character(fData(newobject)[,indexOld])
    
    index <- lapply(newClasses, function(z) 
                    as.vector(sapply(z, function(x) which(featureNames(newobject)==x))))
    
    for (i in 1:length(newClasses)) {
      fData(newobject)[index[[i]], indexNew] <-
        rep(x=names(index)[i], times=length(index[[i]]))
    }
    
    fData(newobject)[,indexNew] <- 
      as.factor(fData(newobject)[,indexNew])
    colnames(fData(newobject))[indexNew] <- newMarkerColumnName
  }
  
  return(newobject)
}

##' Runs the \code{phenoDisco} algorithm.
##' 
##' \code{phenoDisco} is a semi-supervised iterative approach to detect new
##' protein clusters. 
##' 
##'
##' The algorithm performs a phenotype discovery analysis as described in
##' Breckels et al. Using this approach one can identify putative subcellular 
##' groupings in organelle proteomics experiments for more comprehensive 
##' validation in an unbiased fashion. The method is based on the work of Yin 
##' et al. and used iterated rounds of Gaussian Mixture Modelling using the 
##' Expectation Maximisation algorithm combined with a non-parametric outlier
##' detection test to identify new phenotype clusters.
##'
##' One requires 2 or more classes to be labelled in the data 
##' and at a very minimum of 6 markers per class to run the algorithm.
##' The function will check and remove feature with missing values using
##' the \code{\link{filterNA}} method. 
##'
##' Important: Prior to version 1.1.2 the row order in the output was different from
##' the row order in the input. This has now been fixed and row ordering is now
##' the same in both input and output objects.
##' 
##' @param object An instance of class \code{MSnSet}.
##' @param fcol A \code{character} indicating the organellar markers
##' column name in feature meta-data. Default is \code{markers}.
##' @param times Number of runs of tracking. Default is 100.
##' @param GS Group size, i.e how many proteins make a group. Default
##' is 10 (the minimum group size is 4).
##' @param allIter \code{logical}, defining if predictions for all
##' iterations should be saved. Default is \code{FALSE}.
##' @param p Significance level for outlier detection. Default is
##' 0.05.
##' @param seed An optional \code{numeric} of length 1 specifing the
##' random number generator seed to be used.
##' @param verbose Logical, indicating if messages are to be
##' printed out during execution of the algorithm.
##' @return An instance of class \code{MSnSet} containg the \code{phenoDisco}
##' predictions.
##' @author Lisa M. Breckels <lms79@@cam.ac.uk>
##' @references
##' Yin Z, Zhou X, Bakal C, Li F, Sun Y, Perrimon N, Wong ST. Using
##' iterative cluster merging with improved gap statistics to perform online
##' phenotype discovery in the context of high-throughput RNAi screens. BMC
##' Bioinformatics. 2008 Jun 5;9:264.
##' PubMed PMID: 18534020; PubMed Central PMCID: PMC2443381.
##' 
##' Breckels LM, Gatto L, Christoforou A, Groen AJ, Lilley KS and Trotter MWB.
##' The Effect of Organelle Discovery upon Sub-Cellular Protein Localisation. 
##' J Proteomics. In Press.  
##' @examples
##' \dontrun{
##' library(pRolocdata)
##' data(tan2009r1)
##' pdres <- phenoDisco(tan2009r1, fcol = "PLSDA")
##' getPredictions(pdres, fcol = "pd", scol = NULL)
##' plot2D(pdres, fcol = "pd")
##' }
phenoDisco <- function(object,
                       fcol = "markers",
                       times = 100,
                       GS = 10,
                       allIter = FALSE,
                       p = 0.05,
                       seed,
                       verbose = TRUE) {
  ## phenoDisco.R (Lisa's phenoDisco code - last updated 27/02/2013)
  ## fcol = feature column
  ## times = number of runs of tracking 
  ## GS = group size 
  ## p = significance level for outlier detection test

  ## Check GS not outside limits
  if (GS < 4) 
    stop("Group size specified too small")

  if (!missing(seed)) {
    seed <- as.integer(seed)
    set.seed(seed)
  }

  if (!fcol %in% fvarLabels(object))
    stop("'", fcol, "' not found in feature variables.")

  ## Filter data (check no missing values/remove rows with NA)
  if (any(is.na(exprs(object)))) {
      warning("Removing features with missing values.")
      object <- filterNA(object, pNA = 0)
  }
    
  ## Note row order
  fnames <- featureNames(object)
  
  ## Initial settings
  indexFcol <- which(colnames(fData(object)) == fcol)
  test <- table(fData(object)[,indexFcol])
  
  if (any(sapply(test, function(x) x<6)))
    stop("Not enough markers to run phenoDisco: Require > 6 markers per classes")
  
  if (length(test) < 3)
    stop("Not enough classes specified to run phenoDisco: 
          Require a minimum of 2 labelled classes")
  track <- phenotypes <- vector("list")
  currentClasses <- list()
  cond1 <- cond2 <- TRUE
  original <- fcol
  i <- 0
  ## Remove duplicated rows (i.e. identical profiles and add back later)
  duplicatedRows <- FALSE
  if (anyDuplicated(exprs(object))>0) {
    duplicatedRows <- TRUE
    foo <- duplicated(exprs(object))
    duplicateSet <- object[foo,]
    fData(duplicateSet)$pd <- as.character(fData(duplicateSet)$markers)
    uniqueSet <- object[!foo,]
    object <- uniqueSet
  }  
  ## Initiate while loop until no new merges or phenotypes
  while (cond1 & cond2) {
    i <- i+1
    ## ===> Call "tracking.R" to get - 
    ## (1) clusterIDs (from rounds of clustering using GMMs - could use hierarchical)
    ## (2) new members of known classes (from outlier detection using GMMs)
    if (verbose)
      message(paste("Iteration", i)) 
    track[[i]] <- replicate(n = times, 
                            expr = tracking(
                            data = object, 
                            markerCol = fcol,
                            alpha = p))
    ## Update known classes with members assigned to that class
    ## over all iterations of tracking
    classes <- track[[i]][,1]$k
    indK <- list()
    for (j in 1:length(classes)) {
      indK[[j]] <- apply(track[[i]], 2, function(z) which(z$k==classes[j]))
    }
    names(indK) <- classes
    update <- list()
    for (k in 1:length(indK)) {
      tmp <- NULL
      for (j in 1:ncol(track[[i]])) {
        tmp <- c(tmp, rownames(track[[i]][,j]$L[[indK[[k]][j]]]))
      } 
      update[[k]] <- names(which(table(tmp)==ncol(track[[i]])))
    }
    names(update) <- classes
    currentClasses[[i]] <- update
    
    ## Any new phenotypes?
    phenotypes[[i]] <- getNewClusters(track[[i]], groupSize = GS, jc = 1)
    newPhenoName = paste(".pd", i, sep="")
    
    ## CONDITIONS TO STOP LOOP
    if (i != 1) {
      order1 <- order(names(currentClasses[[i-1]]))
      order2 <- order(names(currentClasses[[i]]))
      cond1 <- length(unlist(currentClasses[[i-1]][order1])) !=
        length(unlist(currentClasses[[i]][order2]))
      cond2 <- length(phenotypes[[1]]$coords) > 0
    }
    ## Update MSnsetObject to include new phenotypes as classes
    object <- updateobject(object, 
                           phenotypes[[i]],
                           currentClasses[[i]],
                           oldMarkerColumnName = fcol,
                           newMarkerColumnName = newPhenoName,
                           originalMarkerColumnName = original)
    fcol <- newPhenoName 
  } ## end of while
  foo <- length(names(fData(object)))
  names(fData(object))[foo] <- "pd"
  
  # Add back in any duplicated rows with localisation assigned from pd
  if (duplicatedRows) {
    ind <- apply(exprs(duplicateSet), 1, function(x) 
      which(apply(exprs(object), 1, function(z) all(x==z))))
    
    for (i in 1:length(ind)){
      fData(duplicateSet)$pd[i] <- as.character(fData(object)$pd[ind[i]])
    }
    
    fData(object)$pd <- as.character(fData(object)$pd)
    object <- combine(object, duplicateSet)
  }
  
  if (missing(seed)) {
    procmsg <- paste0("Run phenoDisco using '", original, "': ", date())    
  } else {
    procmsg <- paste0("Run phenoDisco using '", original,
                      "' (seed, ", seed, "): ", date())
  }
  procmsg <- c(procmsg,
               paste0("  with parameters times=", times,
                      ", GS=", GS,
                      ", p=", p, "."))
                      ## ", r=", r, "."))
  object@processingData@processing <-
    c(processingData(object)@processing,
      procmsg)

  if (!allIter) {
      ## FIXME there could be an issue here 
      ## if there were other matching columns
      idx <- grep(".pd", fvarLabels(object))
      fData(object) <- fData(object)[, -idx]
  }
  a <- match(fnames, featureNames(object))
  object <- object[a,]
  object <- MSnbase:::nologging(object, n = 1)
  stopifnot(featureNames(object)==fnames)
  if (validObject(object))
      return(object)
}
