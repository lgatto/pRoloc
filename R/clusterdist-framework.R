######################################################################
## ClustDist: this class summarises the algorithm information from
## running the clustDist algorithm. Information such as (i) the number 
## of k's tested for the kmeans, (ii) the mean and (iii) normalised, 
## pairwise Euclidean distances and (iv) cluster size per numer of 
## component clusters tested, for each GO term tested.
setClass("ClustDist", 
         representation(k = "numeric",
                        dist = "list",
                        term = "character",
                        id = "character",
                        nrow = "numeric",
                        clustsz = "list",
                        components = "vector",
                        fcol = "character"))

setMethod("show",
          signature(object = "ClustDist"),
          function(object) {
            cat("Object of class \"",class(object),"\"\n",sep="")
            cat("fcol = ", object@fcol, "\n")
            cat(" term = ", object@term, "\n")
            cat(" id = ", object@id, "\n")
            cat(" nrow = ", object@nrow, "\n")
            cat("k's tested:", object@k, "\n")
            for (i in 1:length(object@k))
              cat("  Size: ", 
                  paste(object@clustsz[[i]], collapse = ", ")
                  , "\n")
            dfr <- getClusterInfo(object)
            cat("Clusters info:\n")
            print(dfr)
            invisible(NULL)
          })




###################################################################
## ClustDistList: container for multiple ClustDist objects.
.ClustDistList <-
  setClass("ClustDistList",
           slots = c(x = "list",
                     log = "list"),
           contains = "Versioned",
           prototype = prototype(
             new("Versioned",
                 versions = c(ClustDistList = "0.1.0"))),
           validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (!listOf(object@x, "ClustDist"))
               msg <- validMsg(msg, "Not all items are ClustDist objects")
             nvals <- sapply(object@x, validObject)
             if (!all(nvals))
               msg <- validMsg(msg,
                               paste(sum(!nvals),
                                     "ClustDists are not valid."))
             if (is.null(msg)) TRUE
             else msg
           })

ClustDistList <-
  function(x, log = list(call = match.call()))
    .ClustDistList(x = x, log = log)

setMethod("show", "ClustDistList",
          function(object) {
            cat("Instance of class '", class(object), "' containig ",
                length(object), " objects.\n", sep = "")
          })

setMethod("length", "ClustDistList", function(x) length(x@x))

setMethod("names", "ClustDistList", function(x) names(x@x))

setMethod("[", c("ClustDistList", "ANY", "missing", "missing"),
          function(x, i, j = "missing", drop = "missing")
            .ClustDistList(x = msnsets(x)[i]))

setMethod("[[", c("ClustDistList", "ANY", "missing"),
          function(x, i, j = "missing", drop = "missing") {
            if (length(i) != 1)
              stop("subscript out of bounds")
            msnsets(x)[[i]]
          })

setReplaceMethod("names", "ClustDistList",
                 function(x, value) {
                   names(x@x) <- value
                   x
                 })

clustdists <- function(object) object@x

setMethod("lapply", "ClustDistList",
          function(X, FUN, ...) {
            ans <- lapply(clustdists(X), FUN, ...)
            if (listOf(ans, "ClustDist"))
              ans <- ClustDistList(ans)
            ans
          })

setMethod("sapply", "ClustDistList",
          function(X, FUN, ...) {
            ans <- sapply(clustdists(X), FUN, ...)
            if (listOf(ans, "ClustDist"))
              ans <- ClustDistList(ans)
            ans
          })


###################################################################



##      x =  Object of "ClustDistList"
## method =  One of "norm" or "mean", the default is "norm", indicating whether
##           the mean distance, or mean normalised distance per k clusters 
##           should be calculated
##      p =  Normalisation factor. Default is 1/3.
##  nchar = Maximum number of characters of GO ID, before their truncation. Default is 40.
##   ...  = Arguments passed to axis
setMethod("plot", c("ClustDistList", "missing"),
          function(x, y,
                   method = "norm",
                   p = 1/3,
                   main = "",
                   mai = c(5, 1.1, .5, .56),
                   mar = c(15, 6, 2, 2),
                   nchar = 40,
                   ...
          ) {
            opar <- par(no.readonly = TRUE)
            on.exit(par(opar))
            if (!(method == "norm" | method == "mean"))
              stop("method must be one of 'norm' or 'mean' see ?ClustDistList for details")
            if (method == "norm")
              clusterdists <- lapply(x, function(z) normDist(z, best = FALSE, p = p))
            if (method == "mean")
              clusterdists <- lapply(x, function(z) meanDist(z, best = FALSE))
            orgnames <- sapply(x, function(z) z@id)
            ll <- sapply(orgnames, nchar)
            to.shorten <- which(ll > nchar)
            if (length(to.shorten) > 0) {
              replace.with <- paste0(substr(orgnames[to.shorten], 1, 37), "...")
              orgnames[to.shorten] <- replace.with
            }
            min.dist <- sapply(clusterdists, min, na.rm = TRUE)
            oo <- order(min.dist, decreasing = TRUE)
            names(clusterdists) <- orgnames
            orgnames <- orgnames[oo]
            clusterdists <- clusterdists[oo]
            if (method == "norm") 
              ylab = "Mean normalised \nEuclidean distance (per k)"
            else
              ylab = "Mean Euclidean \ndistance (per k)"
            par(mai = c(5, 1.1, .5, .56), mar = c(15, 6, 2, 2) + 0.1)
            boxplot(clusterdists, type = "p", pch = 19, cex = .5, 
                    xaxt = "n", yaxt = "n", axt = "n", xlab = NA, 
                    ylab = ylab, main = main, cex.main = .8)
            axis(1, at = 1:length(oo), labels = orgnames, 
                 las = 2, cex.axis = .7, ...)
            axis(2, at = seq(0, max(unlist(clusterdists), na.rm = TRUE), by = .05), 
                 las = 2, cex.axis = .9)
          })

##      x =  Object of "ClustDist"
##      y =  Object of "MSnSet"
## method =  One of "norm" or "mean", the default is "norm", indicating whether
##           the mean distance, or mean normalised distance per k clusters 
##           should be calculated
##      p =  Normalisation factor. Default is 1/3.
##      nchar = Maximum number of characters of GO ID, before their truncation. Default is 40.
setMethod("plot", c("ClustDist", "MSnSet"),
          function(x, y,
                   method = "norm",
                   p = 1/3,
                   nchar = 40
          ) {
            opar <- par(no.readonly = TRUE)
            on.exit(par(opar))
            if (!(method == "norm" | method == "mean"))
              stop("method must be one of 'norm' or 'mean' see ?plotClustDist for details")
            if (method == "norm")
              clusterdists <- normDist(x, best = FALSE, p = p)
            if (method == "mean")
              clusterdists <- meanDist(x, best = FALSE, p = p)
            clusterdists <- signif(clusterdists, 3)
            numk <- length(x@k)
            if (numk > 3)
              par(mfrow = c(floor(sqrt(numk)), ceiling(sqrt(numk))),
                  oma = c(0, 0, 2, 0))
            else
              par(mfrow = c(1, numk), oma = c(0, 0, 2, 0))
            title <- x@id
            for (i in 1:length(x@k)) {
              fData(y)$.tmp <- "unknown"
              components <- x@components[[i]]
              fData(y)[names(components), ".tmp"] <- as.character(components)
              plot2D(y, ".tmp", main = paste("distance = ", clusterdists[i]))
            }
            mtext(title, outer = TRUE, cex = 1, font = 2)
          })

## Function to output the Euclidean distance matrix and component IDs 
## of each protein for a given k tested
.clusterDistK <- function(k, x, min.size) {
  if (k < nrow(x)) {
    kmcl = kmeans(x, centers = k)
    comps = kmcl$cluster
    ids = tapply(comps, comps, names)
    ll = sapply(ids, length)
    torm = names(which(ll < min.size)) ## Remove components where num of prots < min.size
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

## Calculates Euclidean distances for all columns in a markers matrix 
## input is a `MSnSet` containing the markers matrix. 
## Output is a `EucDist` object with all distances.
clustDist <- function(object,
                     k = 1:5,
                     fcol = "GOMarkers",
                     n = 5,
                     verbose = TRUE) {
  
  ## check min cluster size is not > available GO marker sets
  min.cs <- min(colSums(fData(object)[, fcol]))
  if (min.cs < n)
    stop("There are some columns in fcol = ", fcol, " that have < n proteins.
  Please run filterMinMarkers with n = ", n, " or decrease the size of n.")
  if (min.cs < 2)
      stop("Please run filterMinMarkers. There are some columns in 
       fcol = ", fcol, " with only 1 protein. Can not create a
       cluster with only 1 protein.")
  
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
    cs <- vector("list", length(k))   ## cluster size
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
    termnames = colnames(.pm)
    if (length(grep("GO:", termnames)) > 0) {
      term = termnames
      id = goIdToTerm(term, names = FALSE)
      if (length(grep("/ ", termnames)) > 0) {
        id <- sapply(termnames, function(z) 
          paste(goIdToTerm(strsplit(z, split = "/ ")[[1]], names = FALSE), 
                collapse = "/ "))
        names(id) <- NULL
      }
    } else {
      term = goTermToId(termnames, names = FALSE)
      id = termnames
      if (length(grep("/ ", termnames)) > 0) {
        term <- sapply(termnames, function(z) 
          paste(goTermToId(strsplit(z, split = "/ ")[[1]], names = FALSE), 
                collapse = "/ "))
        names(term) <- NULL
      }
    }
    
    ## Store results in a ClustDist object
    res[[i]] = new("ClustDist", 
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
    dn = dn[which.min(dn)]
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
    dn = dn[which.min(dn)]
  if (k)
    dn <- list(res = dn, k = cc)
  return(dn)
}

## Return numeric of best distances per term
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