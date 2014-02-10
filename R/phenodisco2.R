## tracking.R (L Breckels)
##
## Changes: 2014-02-03 added ndims support
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
tracking <- function(data, alpha = 0.05, markerCol = "markers", ndims=2) {
  ## ===STAGE 1=== TRANSFORM DATA by PCA to reduce data complexity
  k <- names(table(fData(data)[, markerCol]))
  k <- k[which(k!="unknown")]
  k <- sample(k) ## Sample k to avoid bias (order of classes *does* affect cluster 
                 ## members but should not affect the ID of new phenotypes)
  pca <- prcomp(exprs(data), scale=T)$x
  pca <- pca[, 1:ndims]
  tmp <- lapply(k, function(z) fData(data)[, markerCol]==z)
  L <- lapply(tmp, function(z) pca[z , ])
  X <- pca[fData(data)[, markerCol]=="unknown", ]
  originalL <- L
  originalX <- X
  names(originalL) <- names(L) <- k

  ## Define data structure
  track <- newClassData <- vector("list", length(k))
  for (i in 1:length(k)) {
    if (dim(X)[1] > 0) {
      ## ===STAGE 2=== Cluster the set L(k) U X using a GMM and -
      ## (A)  Identify potential new members/candidates of organelle class k
      ## (B)  Get cluster IDs for defining phenotypes at a later stage  
      kUx <- rbind(L[[i]], X)
      gmmCluster <- Mclust(kUx, modelNames = c("EEE","EEV","VEV","VVV"))
      track[[i]] <- gmmCluster$classification # Get cluster number for each prot
      Nclass <- nrow(L[[i]])
      classifyL <- track[[i]][1:Nclass]  
      classifyX <- track[[i]][(Nclass+1):length(track[[i]])] 
      ## Do not include clusters where there is only 1 member! Can not call
      ## outlier detect when cluster contains only one profile! Set min group = 5
      clusterID <- names(which(table(classifyL) >= 5))
      percentL <- sapply(clusterID, function(z) 
        (sum(classifyL == z)/sum(track[[i]] == z)*100))
      ## If one of the clusters contains ONLY labelled members amd no unassigned  
      ## proteins we remove the cluster as outlier detection can not be called!
      if (sum(percentL==100) > 0) {clusterID <- names(which(percentL!=100))}
      ## An if/else statement is here for the instance when we have no
      ## profiles to consider i.e. we get 100% labelled profiles in *ALL* clusters
      if (length(clusterID) < 1) {
        newClassData[[i]] <- L[[i]]; X <- X
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
        ## ===Stage 3=== Call outlier detection on candidates
        outliers <- gmmOutlier(L[[i]], candidates, p=alpha)
        ## ===Stage 4=== Merge/reject candidates. Profiles classified merged
        ## with the current class and rejected members are returned to X
        classified <- names(which(outliers==FALSE))
        if (length(classified) > 0) {
          toAdd <- X[classified, ]
          if(!is.matrix(toAdd)) {
            toAdd <- t(as.matrix(toAdd))
            rownames(toAdd) <- classified
          }
          L[[i]] <- rbind(L[[i]], toAdd)
          torm <- sapply(classified, function(z) which(rownames(X)==z))
          X <- X[-torm, ]
        } else {
          L[[i]] <- L[[i]]
        }
      }
    } else {
      gmmCluster <- Mclust(L[[i]], modelNames = c("EEE","EEV","VEV","VVV"))
      track[[i]] <- gmmCluster$classification
      L[[i]] <- L[[i]]
    }    
  }
  
  list(trackHistory = track, 
       X = X,
       originalX = originalX,
       L = L,
       originalL = originalL,
       k = k)
}


##' Runs the \code{phenoDisco} algorithm.
##' 
##' \code{phenoDisco} is a semi-supervised iterative approach to detect new
##' protein clusters. 
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
##' The function will check and remove features with missing values using
##' the \code{\link{filterNA}} method.
##'
##' A parallel implementation, relying on the \code{BiocParallel}
##' package, has been added in version 1.3.9. See the \code{BPPARAM}
##' arguent for details.
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
##' @param ndims Number of principal components to use as input for
##' the disocvery analysis. Default is 2. Added in version 1.3.9.
##' @param BPPARAM Support for parallel processing using the
##' \code{BiocParallel} infrastructure. When missing (default), the
##' default registered \code{BiocParallelParam} parameters are
##' used. Alternatively, one can pass a valid \code{BiocParallelParam}
##' parameter instance: \code{SnowParam}, \code{MulticoreParam},
##' \code{DoparParam}, \ldots see the \code{BiocParallel} package for
##' details. To revert to the origianl serial implementation, use
##' \code{NULL}.
##' @param seed An optional \code{numeric} of length 1 specifing the
##' random number generator seed to be used. Only relevant when
##' executed in serialised mode with \code{BPPARAM = NULL}. See
##' \code{BPPARAM} for details.
##' @param verbose Logical, indicating if messages are to be
##' printed out during execution of the algorithm.
##' @return An instance of class \code{MSnSet} containing the \code{phenoDisco}
##' predictions.
##' @author Lisa M. Breckels <lms79@@cam.ac.uk>
##' @references
##' Yin Z, Zhou X, Bakal C, Li F, Sun Y, Perrimon N, Wong ST. Using
##' iterative cluster merging with improved gap statistics to perform
##' online phenotype discovery in the context of high-throughput RNAi
##' screens. BMC Bioinformatics. 2008 Jun 5;9:264.  PubMed PMID:
##' 18534020.
##' 
##' Breckels LM, Gatto L, Christoforou A, Groen AJ, Lilley KS and
##' Trotter MWB.  The Effect of Organelle Discovery upon Sub-Cellular
##' Protein Localisation.  J Proteomics. 2013 Aug 2;88:129-40. doi:
##' 10.1016/j.jprot.2013.02.019. Epub 2013 Mar 21.  PubMed PMID:
##' 23523639.
##' @examples
##' \dontrun{
##' library(pRolocdata)
##' data(tan2009r1)
##' pdres <- phenoDisco(tan2009r1, fcol = "PLSDA")
##' getPredictions(pdres, fcol = "pd", scol = NULL)
##' plot2D(pdres, fcol = "pd")
##' }
phenoDisco2 <- function(object,
                        fcol = "markers",
                        times = 100,
                        GS = 10,
                        allIter = FALSE,
                        p = 0.05,
                        ndims = 2,
                        BPPARAM,
                        seed,
                        verbose = TRUE) {
    ## phenoDisco.R 
    ## Changes:
    ##   ndims  2014-02-03
    ##   BPARAM 2014-02-03

  ## Check data and parameters properly specified
  if (GS < 4) 
      stop("Group size specified too small.")
  if (!anyUnknown(object))
      stop("No unlabelled features (conventionally marked 'unknown') in your data.")
  if (!missing(seed) & is.null(BPPARAM)) {
    seed <- as.integer(seed)
    set.seed(seed)
  }
  if (!fcol %in% fvarLabels(object))
    stop("'", fcol, "' not found in feature variables.")
  if (any(is.na(exprs(object)))) {
      warning("Removing features with missing values.")
      object <- filterNA(object, pNA = 0)
  }
  ## Check ndims is sensible
  .pca <- prcomp(exprs(object), center = TRUE, scale = TRUE)$x
  if (ndims > ncol(.pca)) {
    warning("ndims > number of principal components available, using maximum 
          number of components (ndims = ", ncol(.pca), ")")
    ndims <- ncol(.pca)
  }
  if (ndims <= 1) {
    warning("ndims <= 1, using ndims = 2")
    ndims <- 2
  }
  ## Check we have enough labelled data to start
  test <- table(fData(object)[,fcol])
  if (any(sapply(test, function(x) x<6)))
    stop("Not enough markers to run phenoDisco: Require > 6 markers per classes")
  if (length(test) < 3)
    stop("Not enough classes specified to run phenoDisco: 
          Require a minimum of 2 labelled classes")
  ## Remove duplicated rows (i.e. identical profiles and add back later)
  duplicatedRows <- FALSE
  if (anyDuplicated(exprs(object))>0) {
    duplicatedRows <- TRUE
    foo <- duplicated(exprs(object))
    duplicateSet <- object[foo,]
    fData(duplicateSet)$pd <- as.character(fData(duplicateSet)[, fcol])
    uniqueSet <- object[!foo,]
    object <- uniqueSet
  }
  
  ## Note row order
  fnames <- featureNames(object)
  
  ## Initial settings
  track <- phenotypes <- vector("list")
  currentClasses <- list()
  cond1 <- cond2 <- TRUE
  original <- fcol
  i <- 0

  ## Initiate while loop until no new merges or phenotypes
  while (cond1 & cond2) {
    i <- i+1
    ## ===> Call "tracking.R" to get - 
    ## (1) clusterIDs (from rounds of clustering using GMMs - could use hierarchical)
    ## (2) new members of known classes (from outlier detection using GMMs)
    if (verbose)
      message(paste("Iteration", i)) 

    if (missing(BPPARAM)) {
        ## default: taking first registered BiocParallelParam
        track[[i]] <- simplify2array(bplapply(seq_len(times),
                                              function(x) 
                                              tracking(data = object,
                                                       alpha = p,
                                                       markerCol = fcol,
                                                       ndims = ndims)))
    } else if (inherits(BPPARAM, "BiocParallelParam")) {
        ## using user-specified BiocParallelParam
        track[[i]] <- simplify2array(bplapply(seq_len(times),
                                              function(x) 
                                              tracking(data = object,
                                                       alpha = p,
                                                       markerCol = fcol,
                                                       ndims = ndims), 
                                                       BPPARAM = BPPARAM))
    } else if (is.null(BPPARAM)) {
        ## serialised version (original implementation)
        track[[i]] <- replicate(n = times, 
                                expr = tracking(
                                    data = object, 
                                    markerCol = fcol,
                                    alpha = p,
                                    ndims = ndims))
    } else {
        stop("Non valid BPPARAM. See ?phenoDisco for details.")
    }
    
    ## Update known classes with members assigned to that class
    ## over all iterations of tracking
    classes <- track[[i]][,1]$k
    indK <- list()
    for (j in 1:length(classes)) {
      indK[[j]] <- apply(track[[i]], 2, function(z) which(z$k==classes[j]))
    }
    names(indK) <- classes
    update <- vector("list", length(indK))
    for (k in 1:length(indK)) {
      tmp <- NULL
      for (j in 1:ncol(track[[i]])) {
        tmp <- c(tmp, rownames(track[[i]][,j]$L[[indK[[k]][j]]]))
      } 
      update[[k]] <- names(which(table(tmp)==ncol(track[[i]])))
    }
    names(update) <- classes
    currentClasses[[i]] <- update
    
    ## Update X members (consistently X in all replicates)
    idX <- names(which(table(unlist(lapply(track[[i]]["X", ], 
                                           function(z) rownames(z))))==times))
    .tmp <- track[[i]]["originalX", 1][[1]]
    indX <- sapply(idX, function(z) which(rownames(.tmp)==z))
    candidatesX <- .tmp[indX, ]
    
    ## Get track history for X members
    trackingStats <- track[[i]][1,]
    trackingStats <- unlist(trackingStats)
    idHist <- sapply(idX, function(z) trackingStats[which(names(trackingStats) == z)])
    idHist <- t(idHist)
    colnames(idHist) <- NULL
    
    ## Any new phenotypes?
    phenotypes[[i]] <- getNewClusters(idHist, candidatesX, groupSize = GS, jc = 1)
    newPhenoName <- paste(".pd", i, sep="")
    
    ## CONDITIONS TO STOP LOOP
    if (i != 1) {
      order1 <- order(names(currentClasses[[i-1]]))
      order2 <- order(names(currentClasses[[i]]))
      cond1 <- length(unlist(currentClasses[[i-1]][order1])) !=
        length(unlist(currentClasses[[i]][order2]))
      cond2 <- length(phenotypes[[1]]$coords) > 0
    }
    ## Update MSnSetObject to include new phenotypes as classes
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
  procmsg <- paste0(procmsg,
                    "\n   with parameters times=", times,
                    ", GS=", GS,
                    ", p=", p,
                    ", ndims=", ndims, ".")    
  object <- MSnbase:::logging(object,  procmsg, date. = FALSE)

  if (!allIter) {
      ## FIXME there could be an issue here 
      ## if there were other matching columns
      idx <- grep(".pd", fvarLabels(object))
      fData(object) <- fData(object)[, -idx]
  }
  a <- match(fnames, featureNames(object))
  object <- object[a,]
  object <- MSnbase:::nologging(object, n = 1)
  stopifnot(featureNames(object) == fnames)
  if (validObject(object))
      return(object)
}
