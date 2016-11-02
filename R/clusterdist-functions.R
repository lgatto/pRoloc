
## Sub-function for clustDist: output the Euclidean distance matrix and 
## component IDs of each protein for a given k tested
.clusterDistK <- function(k, x, min.size) {
  if (k < nrow(x)) {
    kmcl <- kmeans(x, centers = k)
    comps <- kmcl$cluster
    ids <- tapply(comps, comps, names)
    ll <- sapply(ids, length)
    torm <- names(which(ll < min.size)) ## Remove components where num of prots < min.size
    if (length(torm) > 0) {
      if (length(ll) != length(torm)) {
        indrm <- sapply(torm, function(x) which(names(ids) == x))
        ids <- ids[-indrm]
      }
    }
    if (length(ll) != length(torm)) {
      res <- lapply(ids,
                    function(z) dist(x[z, ]))
      list(res = res,
           comps = comps)
    } else {
      list(res = NA,
           comps = NA)
    }
  } else {
    list(res = NA,
         comps = NA)
  }
}


##' Pairwise Distance Computation for Protein Information Sets 
##' 
##' This function computes the mean (normalised) pairwise
##' distances for pre-defined sets of proteins. 
##' 
##' The input to the function is a \code{MSnSet} dataset 
##' containing a matrix appended to the feature data slot 
##' identifying the membership of protein instances to 
##' a pre-defined set(s) e.g. a specific Gene Ontology term etc.
##' 
##' For each protein set, the \code{clustDist} function (i) 
##' extracts all instances belonging to the set, (ii) using 
##' the \code{kmeans} algorithm fits and tests \code{k = c(1:5)} 
##' (default) cluster components to each set, (iii) calculates 
##' the mean pairwise distance for each \code{k} tested.
##' 
##' Note: currently distances are calcualted in Euclidean space, 
##' but other distance metrics will be supported in the future). 
##' 
##' The output is a \code{list} of \code{ClustDist} objects, 
##' one per information cluster. The \code{ClustDist} 
##' class summarises the algorithm information such as the number of k's 
##' tested for the kmeans, and mean and normalised pairwise Euclidean 
##' distances per numer of component clusters tested. See \code{?ClustDist}
##' for more details.
##' 
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param k The number of clusters to try fitting to the protein set. 
##' Default is \code{k = 1:5}.
##' @param fcol The feature meta-data containing matrix of protein sets/
##' marker definitions. Default is \code{GOAnnotations}.
##' @param n The minimum number of proteins per set. If protein sets
##' contain less than \code{n} instances they will be ignored. 
##' Defualt is 5.
##' @param verbose A logical defining whether a progress bar is displayed.
##' @param seed An optional seed for the random number generator.
##' @return An instance of \code{"\linkS4class{ClustDistList}"} containing
##' a \code{"\linkS4class{ClustDist}"} instance for every protein set, which
##' summarises the algorithm information such as the number of k's tested 
##' for the kmeans, and mean and normalised pairwise Euclidean distances 
##' per numer of component clusters tested. 
##' @seealso For class definitions see \code{"\linkS4class{ClustDistList}"} 
##' and \code{"\linkS4class{ClustDist}"}.
##' @author Lisa Breckels
##' @examples
##' library(pRolocdata)
##' data(dunkley2006)
##' par <- setAnnotationParams(inputs =
##'                    c("Arabidopsis thaliana genes",
##'                    "Gene stable ID"))
##' ## add protein sets/annotation information
##' xx <- addGoAnnotations(dunkley2006, par)
##' ## filter
##' xx <- filterMinMarkers(xx, n = 50)
##' xx <- filterMaxMarkers(xx, p = .25)
##' ## get distances for protein sets 
##' dd <- clustDist(xx)
##' ## plot clusters for first 'ClustDist' object 
##' ## in the 'ClustDistList'
##' plot(dd[[1]], xx)
##' ## plot distances for all protein sets 
##' plot(dd)
##' ## Extract normalised distances
##' ## Normalise by n^1/3
##' minDist <- getNormDist(dd, p = 1/3)
##' ## Get new order according to lowest distance
##' o <- order(minDist)
##' ## Re-order GOAnnotations 
##' fData(xx)$GOAnnotations <- fData(xx)$GOAnnotations[, o]
##' if (interactive()) {
##' pRolocVis(xx, fcol = "GOAnnotations")
##' }
clustDist <- function(object,
                      k = 1:5,
                      fcol = "GOAnnotations",
                      n = 5,
                      verbose = TRUE,
                      seed) {
  ## check min cluster size is not > available GO marker sets
  min.cs <- min(colSums(fData(object)[, fcol]))
  if (min.cs < n)
    stop("There are some columns in fcol = ", fcol, " that have < n proteins.
         Please run filterMinMarkers with n = ", n, " or decrease the size of n.")
  if (min.cs < 2)
    stop("Please run filterMinMarkers. There are some columns in 
         fcol = ", fcol, " with only 1 protein. Can not create a
         cluster with only 1 protein.")

  ## allow setting of seed to re-produce results
  if (missing(seed)) {
    seed <- sample(.Machine$integer.max, 1)
  }
  .seed <- as.integer(seed)
  set.seed(.seed)
  
  ## Matrix of potential markers
  .pm <- fData(object)[, fcol]

  ## Progress bar
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = ncol(.pm), style = 3)
    ._k <- 0
  }
  
  ## For each term i in markers matrix
  res <- comps <- vector("list", length = ncol(.pm))
  for (i in 1:ncol(.pm)) {
    
    if (verbose) {
      setTxtProgressBar(pb, ._k)
      ._k <- ._k + 1
    }
    
    ## Get expression data for term i
    data <- exprs(object)[.pm[, i] == 1, ]
    
    ## (i) Try fitting k components, for each k tested calculate Euclidean
    ## distance matrix of each component, (ii) note down number of proteins
    ## per component. (iii) note down membership of each protein in each
    ## component
    dist.res <- vector("list", length(k))
    for (m in 1:length(k)) {
      dist.res[[m]] <- .clusterDistK(k = k[m], data, min.cs)
    }
    eucl <- lapply(dist.res, function(z) z$res)   ## Euclidean dist matrix
    comps <- lapply(dist.res, function(z) z$comps)  ## component ID
    cs <- vector("list", length(k)) ## cluster size
    for (.cs in 1:length(cs)) {
      .tmp <- sapply(eucl[[.cs]], nrow)
      if (is.null(.tmp[[1]]))
        cs[[.cs]] <- NA
      else
        cs[[.cs]] <- .tmp
    }
    ## NB: Any k's that can't be fit e.g. if k's to test are k = 1:6 and the 
    ## minimum cluster size is k = 4  then we will have NA's for k = 5:6
    
    ## GO term names and IDs
    termnames <- colnames(.pm)
    if (length(grep("GO:", termnames)) > 0) {
      term <- termnames
      id <- goIdToTerm(term, names = FALSE, keepNA = FALSE)
      if (length(grep("/ ", termnames)) > 0) {
        id <- sapply(termnames, function(z) 
          paste(goIdToTerm(strsplit(z, split = "/ ")[[1]], names = FALSE, keepNA = FALSE),
                collapse = "/ "))
        names(id) <- NULL
      }
    } else {
      term <- goTermToId(termnames, names = FALSE, keepNA = FALSE)
      id <- termnames
      if (length(grep("/ ", termnames)) > 0) {
        term <- sapply(termnames, function(z)
          paste(goTermToId(strsplit(z, split = "/ ")[[1]], names = FALSE, keepNA = FALSE),
                collapse = "/ "))
        names(term) <- NULL
      }
    }
    
    ## Store results in a ClustDist object
    res[[i]] <- new("ClustDist",
                    k = k, dist = eucl, fcol = fcol, 
                    nrow = nrow(data), clustsz = cs, 
                    components = comps, term = term[i], 
                    id = id[i])
  }
  if (verbose) {
    setTxtProgressBar(pb, ._k)
    close(pb)
  }
  names(res) <- sapply(res, function(x) x@id)
  res <- ClustDistList(res)
  if(validObject(res))
    return(res)
}


## For each k tested get associated distances and take mean distance per k e.g.
## for k = 3, if there are 3 clusters, take Euclidean distance of each cluster
## then calculate the mean of each cluster, then take the mean of all clusters
## for each k
meanDist <- function(object, best = FALSE, k = FALSE) {
  ## calculate mean(di), i = 1..5
  clustinfo <- object@dist
  d <- lapply(1:length(object@k),
              function(z)
                sapply(clustinfo[[z]], mean))
  .indna <- is.na(d)
  cc <- sapply(d, length)
  cc[.indna] <- NA
  dn <- sapply(d, mean)
  if (best)
    dn <- dn[which.min(dn)]
  if (k)
    dn <- list(res = dn, k = cc)
  return(dn)
}



## Calculate normalised distances for each k tested. 
## Let If n = number of proteins per cluster.
## Then we normalise by n^p. Heuristically we find that
## p = 1/3 works best.
normDist <- function(object, ## ClustDist object
                     best = FALSE,
                     k = FALSE,
                     p = 1/3) {
  ## calculate mean(mean(di)/sqrt(ni))
  clustinfo <- object@dist
  d <- lapply(1:length(object@k),
              function(z)
                sapply(clustinfo[[z]], mean))
  cc <- sapply(d, length)
  idx <- sapply(d, names)
  .indna <- is.na(d)
  cc[.indna] <- idx[.indna] <- NA
  ll <- lapply(object@components, table)
  ll <- lapply(1:length(object@k), function(z) ll[[z]][idx[[z]]])
  dn <- lapply(1:length(object@k),
               function(z)
                 d[[z]]/(ll[[z]]^p))
  ## now take mean(mean(di)/ni), i = 1..5
  dn <- sapply(dn, mean)
  if (best)
    dn <- dn[which.min(dn)]
  if (k)
    dn <- list(res = dn, k = cc)
  return(dn)
}


##' Extract Distances from a \code{"ClustDistList"} object
##' 
##' This function computes and outputs normalised distances from a
##' \code{"\linkS4class{ClustDistList}"} object.
##' 
##' @param object An instance of class \code{"\linkS4class{ClustDistList}"}.
##' @param p The normalisation factor. Default is 1/3.
##' @return An numeric of normalised distances, one per protein set in the
##' \code{ClustDistList}.
##' @seealso \code{"\linkS4class{ClustDistList}"}, \code{"\linkS4class{ClustDist}"},
##' and examples in \code{clustDist}.
##' @author Lisa Breckels
getNormDist <- function(object, ## EucDist class
                        p = 1/3) {
  dists <- sapply(object, normDist, best = TRUE, k = FALSE, p = p)
  names(dists) <- sapply(object, function(z) z@id)
  return(dists)
}

## Wrap up results in a data.frame
getClusterInfo <- function(object) {
  ## mean of distances per component per k 
  xx <- meanDist(object, k = TRUE)
  ks.mean <- xx$k
  res.mean <- xx$res
  best <- which.min(res.mean)
  res.mean <- format(res.mean, digits = 4)
  res.mean[best] <- paste0("*", res.mean[best])
  ## normalise by n^1/3
  pp <- normDist(object, k = TRUE, p = 1/3)
  ks.norm <- pp$k
  res.norm <- pp$res
  best <- which.min(res.norm)
  res.norm <- format(res.norm, digits = 4)
  res.norm[best] <- paste0("*", res.norm[best])
  dfr <- data.frame(ks.mean = ks.mean,
                    mean = res.mean,
                    ks.norm = ks.norm,
                    norm = res.norm)
  rownames(dfr) <- paste("k =", object@k)
  return(dfr)
}
