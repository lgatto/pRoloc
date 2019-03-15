##' These functions implement the T augmented Gaussian mixture (TAGM)
##' model with Novelty detection for mass spectrometry-based spatial proteomics datasets.
##' Markov-chain Monte-Carlo (MCMC) is used for inference.
##'
##' The `tagmMcmcTrain_Nov` function generates the samples from the
##' posterior distributions (object of class `MCMCParams`) based on an
##' annotated quantitative spatial proteomics dataset (object of class
##' [`MSnbase::MSnSet`]). Both are then passed to the `tagmNoveltyProcess`
##' function to predict detect new sub-cellular structures. 
##' See the *pRoloc-bayesian-Novelty* vignette for
##' details and examples. In this implementation, if numerical instability
##' is detected in the covariance matrix of the data a small multiple of
##' the identity is added. A message is printed if this conditioning step
##' is performed.
##'
##' @title Localisation of protiens and Novelty detection TAGM MCMC method
##' @param object An [`MSnbase::MSnSet`] containing the spatial
##'     proteomics data to be passed to `tagmMcmcTrain_Nov` and
##'     `tagmNoveltyProcess`.
##' @param fcol The feature meta-data containing marker definitions.
##'     Default is `markers`.
##' @param method A `charachter()` describing the inference method for
##'     the TAGM algorithm. Default is `"MCMC"`.
##' @param numIter The number of iterations of the MCMC
##'     algorithm. Default is 1000.
##' @param burnin The number of samples to be discarded from the
##'     begining of the chain. Default is 100.
##' @param thin The thinning frequency to be applied to the MCMC
##'     chain.  Default is 5.
##' @param mu0 The prior mean. Default is `colMeans` of the expression
##'     data.
##' @param lambda0 The prior shrinkage. Default is 0.01.
##' @param nu0 The prior degreed of freedom. Default is
##'     `ncol(exprs(object)) + 2`
##' @param S0 The prior inverse-wishart scale matrix. Empirical prior
##'     used by default.
##' @param beta0 The prior Dirichlet distribution
##'     concentration. Default is 1 for each class.
##' @param u The prior shape parameter for Beta(u, v). Default is 2
##' @param v The prior shape parameter for Beta(u, v). Default is 10.
##' @param numChains The number of parrallel chains to be run. Default
##'     it 4.
##' @param BPPARAM Support for parallel processing using the
##'     `BiocParallel` infrastructure. When missing (default), the
##'     default registered `BiocParallelParam` parameters are
##'     used. Alternatively, one can pass a valid `BiocParallelParam`
##'     parameter instance: `SnowParam`, `MulticoreParam`,
##'     `DoparParam`, ... see the `BiocParallel` package for
##'     details.
##' @param K_bar The maximum number of new phenotypes allowed to be detected.
##' The default it 5.
##' @return `tagmMcmcTrain` returns an instance of class
##'     `MCMCParams`.
##' @md
##' @references 
##' @rdname tagm-Novelty

tagmMcmcTrain_Nov <- function(object,
                              fcol = "markers",
                              method = "MCMC",
                              numIter = 1000L,
                              burnin = 100L,
                              thin = 5L,
                              mu0 = NULL,
                              lambda0 = 0.01,
                              nu0 = NULL,
                              S0 = NULL,
                              beta0 = NULL,
                              u = 2,
                              v = 10,
                              numChains = 4L,
                              K_bar = 5,
                              BPPARAM = BiocParallel::bpparam()) {
  
  ## get expression marker data
  markersubset <- markerMSnSet(object, fcol = fcol)
  markers <- getMarkerClasses(markersubset, fcol = fcol)
  mydata <- exprs(markersubset)
  X <- exprs(unknownMSnSet(object, fcol = fcol))
  
  ## get data dize
  N <- nrow(mydata)
  D <- ncol(mydata)
  K <- as.integer(length(markers) + K_bar)
  
  ## set priors
  if (is.null(nu0))
    nu0 <- D + 2
  
  if (is.null(S0))
    S0 <- diag( colSums(( mydata - mean( mydata)) ^ 2) / N)/( K ^ (1/D))
  
  if (is.null(mu0))
    mu0 <- colMeans( mydata)
  
  if (is.null(beta0))
    beta0 <- rep(1, K)
  
  ## Store Priors
  .priors <- list(mu0 = mu0,
                  lambda0 = lambda0,
                  nu0 = nu0,
                  S0 = S0,
                  beta0 = beta0)
  
  ## chains run in parallel, repeating number of iterations
  .res <- bplapply(rep(numIter, numChains),
                   FUN = tagmMcmcChain_Nov,
                   object = object,
                   fcol = fcol,
                   method = "MCMC",
                   burnin = burnin,
                   thin = thin,
                   mu0 = mu0,
                   lambda0 = lambda0,
                   nu0 = nu0,
                   S0 = S0,
                   beta0 = beta0,
                   u = u,
                   v = v,
                   K_bar = K_bar,
                   BPPARAM = BPPARAM)
  
  ## Construct class MCMCChains
  .ans <- pRoloc:::.MCMCChains(chains = .res)
  
  ## Construct class MCMCParams
  out <- pRoloc:::.MCMCParams(method = "TAGM.MCMC",
                              chains = .ans,
                              priors = .priors,
                              summary = pRoloc:::.MCMCSummary())
  
  return(out)
}

##' @param object An [`MSnbase::MSnSet`] containing the spatial
##'     proteomics data to passed to `tagmMcmcTrain_Nov`
##' @param params An instance of class `MCMCParams`, as generated by
##'     [tagmMcmcTrain_Nov()]    
##' @param fcol The feature meta-data containing marker definitions.
##'     Default is `markers`.
##' @return `tagmNoveltyChainProcess` returns an instance of class [NoveltyTagm()].

tagmNoveltyProcess <- function(object,
                               params,
                               fcol = "markers") {
  ## Checks that params is object of class MCMCChain
  stopifnot(inherits(params, "MCMCParams"))
  numChains <- length(params)
  
  ## Storage
  .res <- vector("list", length = numChains)
  
  ## Marker Classes
  markersubset <- markerMSnSet(object = object, fcol)
  K_markers <- nlevels(fData(markersubset)[, fcol])
  
  ## For each Chain compute psm and maxpear
  for (j in seq_len(numChains)){
    .res[[j]] <- tagmNoveltyChainProcess(object = object,
                                    params = pRoloc:::chains(params)[[j]],
                                    fcol = fcol)
  }
  
  ## Pool Chains and get psm and maxpear for combined chains
  pooledparams <- mcmc_pool_chains(params)
  .pooledres <- tagmNoveltyChainProcess(object = object,
                                        params = pRoloc:::chains(pooledparams)[[1]],
                                        fcol = fcol)

  # Use pooled params to compute discovery probability
  pooledparams <- tagmMcmcProcess(pooledparams)
  tagm.newcluster.prob <- 1 - rowSums(pooledparams@summary@tagm.joint[, 1:K_markers])
  
  .out <- .NoveltyChains(pooledRes = .pooledres,
                         NoveltyChains = .res,
                         tagm.newcluster.prob = tagm.newcluster.prob)
  
  
  return(.out)
}

##' @param object An [`MSnbase::MSnSet`] containing the spatial
##'     proteomics data to be passed to `tagmMcmcTrain_Nov` and
##'     `tagmNoveltyProcess`.
##' @param params An instance of class `MCMCParams`, as generated by
##'     [tagmMcmcTrain_Nov()]    
##' @param fcol The feature meta-data containing marker definitions.
##'     Default is `markers`.
##' @return `tagmNoveltyChainProcess` returns an instance of class [NoveltyChain()].`

tagmNoveltyChainProcess <- function(object,
                                    params,
                                    fcol = "markers") {
  ## Checks that params is object of class MCMCChain
  stopifnot(inherits(params, "MCMCChain"))
  
  ## Make marker allocation matrix and align rownames
  markersubset <- markerMSnSet(object = object, fcol)
  markers <- getMarkerClasses(markersubset, fcol = fcol)
  markerallocationMat <- matrix(NA, ncol = ncol(params@Component), nrow = nrow(markersubset))
  rownames(markerallocationMat) <- rownames(markersubset)
  
  K <- params@K
  K_markers <- nlevels(fData(markersubset)[, fcol])
  
  
  ## populate marker allocation matrix
  for (j in seq.int(K_markers)) {
    toSubset <- rownames(markersubset)[fData(markersubset)[, fcol] == markers[j]]
    markerallocationMat[toSubset, ] <- rep(j, ncol(params@Component))
  }
  
  # Form component matrix with markers
  combinedComponent <- rbind(params@Component, markerallocationMat)
  
  #Compute PSM and summarise using MCclust methods
  modal.k <- which.max(tabulate(apply(params@Component, 2, function(x) length(unique(x)))))
  psm <- mcclust::comp.psm(t(combinedComponent))
  mp <- mcclust::maxpear(psm = psm, max.k = K)
  
  # Match rownames 
  names(mp$cl) <- rownames(combinedComponent)
  rownames(psm) <- colnames(psm) <- c(rownames(unknownMSnSet(object)), rownames(markersubset))
  
  # Finding mapping from cluster to organelle because of label switching
  mapping <- matrix(NA, ncol = K_markers, nrow = 1)
  for (j in seq.int(K_markers)) {
    toSubset <- rownames(combinedComponent
                         [rownames(markersubset), ])[combinedComponent[rownames(markersubset), 1] == j]
    mapping[j] <- mp$cl[toSubset][1]
  }
  mapping <- c(mapping, seq.int(1:K)[!seq.int(1:K) %in% mapping])
  
  ## Match up clusters with organelles and phenotypes
  for (j in seq.int(K)) {
    mp$cl[mp$cl == mapping[j]] <- rownames(params@ComponentParam@mk)[j]
  }
  
  ## Make factor
  mp$cl <- factor(mp$cl)
  
  ## Create maxPear class
  maxPear <- .maxPear(cl = mp$cl,
                      value = mp$value,
                      method = mp$method)
  
  ## Make novelty Chain
  .out <- .NoveltyChain(psm = psm,
                         mp = maxPear)
  
  return(.out) 
}


##' @param object An [`MSnbase::MSnSet`] containing the spatial
##'     proteomics data to be passed to `tagmMcmcTrain` and
##'     `tagmPredict`.
##' @param fcol The feature meta-data containing marker definitions.
##'     Default is `markers`.
##' @param method A `charachter()` describing the inference method for
##'     the TAGM algorithm. Default is `"MCMC"`.
##' @param numIter The number of iterations of the MCMC
##'     algorithm. Default is 1000.
##' @param burnin The number of samples to be discarded from the
##'     begining of the chain. Default is 100.
##' @param thin The thinning frequency to be applied to the MCMC
##'     chain.  Default is 5.
##' @param mu0 The prior mean. Default is `colMeans` of the expression
##'     data.
##' @param lambda0 The prior shrinkage. Default is 0.01.
##' @param nu0 The prior degreed of freedom. Default is
##'     `ncol(exprs(object)) + 2`
##' @param S0 The prior inverse-wishart scale matrix. Empirical prior
##'     used by default.
##' @param beta0 The prior Dirichlet distribution
##'     concentration. Default is 1 for each class.
##' @param u The prior shape parameter for Beta(u, v). Default is 2
##' @param v The prior shape parameter for Beta(u, v). Default is 10.
##' @param K_bar The maximum number of new phenotypes allowed to be detected.
##' The default it 5.
##' @return `tagmMcmcChain` returns an instance of class
##'     [MCMCChain()].
##' @md
##' @noRd
tagmMcmcChain_Nov <- function(object,
                              fcol = "markers",
                              method = "MCMC",
                              numIter = 1000L,
                              burnin = 100L,
                              thin = 5L,
                              mu0 = NULL,
                              lambda0 = 0.01,
                              nu0 = NULL,
                              S0 = NULL,
                              beta0 = NULL,
                              u = 2,
                              v = 10,
                              K_bar = 5) {
  if (burnin >= numIter)
    stop("Burnin is larger than iterations you will not retain any samples")
  
  ## on the fly number of samples to be retained
  retained <- seq.int(burnin + 1L, numIter , thin)
  
  ## get expression marker data
  markersubset <- markerMSnSet(object, fcol = fcol)
  markers <- getMarkerClasses(markersubset, fcol = fcol)
  mydata <- exprs(markersubset)
  X <- exprs(unknownMSnSet(object, fcol = fcol))
  
  ## get data dize
  N <- nrow(mydata)
  D <- ncol(mydata)
  K_markers <- length(markers)
  K <- as.integer(length(markers) + K_bar)
  
  ## set empirical priors
  if (is.null(nu0))
    nu0 <- D + 2
  
  if (is.null(S0))
    S0 <- diag( colSums(( mydata - mean( mydata)) ^ 2) / N)/( K ^ (1/D))
  
  if (is.null(mu0))
    mu0 <- colMeans( mydata)
  
  if (is.null(beta0))
    beta0 <- rep(1, K)
  
  ## save priors
  .priors <- list(mu0 = mu0,
                  lambda0 = lambda0,
                  nu0 = nu0,
                  S0 = S0,
                  beta0 = beta0)
  
  ## create storage for posterior parameters
  mk <- matrix(0, nrow = K, ncol = D)
  lambdak <- matrix(0, nrow = K, ncol = 1)
  nuk <- matrix(0, nrow = K, ncol = 1)
  sk <- array(0, dim = c(K, D, D))
  
  ## create storage for cluster parameters
  xk <- matrix(0, nrow = K, ncol = D)
  
  ## update prior with training data
  nk <- tabulate(fData(markersubset)[, fcol])
  nk <- c(nk, rep(0, K_bar)) # attach pheno-discovery clusters
  
  for (j in seq.int(K_markers))
    xk[j, ] <- colSums(mydata[fData(markersubset)[, fcol] == markers[j], ])/nk[j]
  
  lambdak <- lambda0 + nk
  nuk <- nu0 + nk
  mk <- t((t(nk * xk) + lambda0 * mu0)) / lambdak
  betak <- beta0 + nk
  
  for(j in seq.int(K)) {
    if (j > K_markers) { ## pheno-discovery clusters get prior
      sk[j, , ] <- S0
    } else {
      sk[j, , ] <- S0 + t(mydata[fData(markersubset)[, fcol] == markers[j], ]) %*%
        mydata[fData(markersubset)[, fcol] == markers[j],] +
        lambda0 * mu0 %*% t(mu0) - lambdak[j] * mk[j, ] %*% t(mk[j, ])
    }
  }
  
  ## global parameters
  M <- colMeans(exprs(object))
  V <- cov(exprs(object))/2
  eigsV <- eigen(V)
  if (min(eigsV$values) < .Machine$double.eps) {
    V <- cov(exprs(object))/2 + diag(10^{-6}, D)
    message("co-linearity detected; a small multiple of
            the identity was added to the covariance")
  }
  
  ## storage
  Component <- matrix(0, nrow = nrow(X), ncol = numIter)
  ComponentProb <- array(0, c(nrow(X), numIter, K))
  Outlier <- matrix(0, nrow = nrow(X), ncol = numIter)
  OutlierProb <- array(0, c(nrow(X), numIter, 2))
  
  ## initially assigned all unlabelled points to clusters greedily
  for(j in seq.int(K))
    ComponentProb[, 1, j] <- pRoloc:::dmvtCpp(X,
                                              mu_ = mk[j, ],
                                              sigma_ = (1 + lambdak[j]) * sk[j, , ] / (lambdak[j] * (nuk[j] - D + 1)),
                                              df_ = nuk[j] - D + 1,
                                              log_ = TRUE,
                                              ncores_ = 1,
                                              isChol_ = FALSE)
  
  Component[, 1] <- apply(X = ComponentProb[, 1, ], 1, FUN = which.max)
  
  ## randomly reassign half the proteins to empty clusters
  prm <- sample(x = seq.int(1:nrow(X)), size = floor(nrow(X)/2), replace = FALSE)
  Component[prm, 1] <- sample(x = seq.int(1:K_bar), size = floor(nrow(X)/2), replace = TRUE) + K_markers
  
  ## Need statistics for pheno clusters
  nk <- nk + tabulate(Component[prm, 1])
  
  for (j in seq.int(K_markers + 1, K))
    xk[j, ] <- colSums(X[prm,][Component[prm, 1] == j,])/nk[j]
  
  lambdak <- lambda0 + nk
  nuk <- nu0 + nk
  mk <- t((t(nk * xk) + lambda0 * mu0)) / lambdak
  betak <- beta0 + nk
  
  for (j in seq.int(K_markers + 1, K))
    sk[j, , ] <- S0 + t(X[prm,][Component[prm, 1] == j,]) %*% X[prm,][Component[prm, 1] == j,] +
    lambda0 * mu0 %*% t(mu0) - lambdak[j] * mk[j, ] %*% t(mk[j, ])
  
  
  ## initially assign all proteins to global component unless they've been assigned to new phenotype
  Outlier[seq.int(nrow(X))[-prm], 1] <- 0
  Outlier[prm, 1] <- 1
  
  ## initial allocation statistics
  tau1 <- sum(Outlier[, 1] == 1) + N
  tau2 <- sum(Outlier[, 1] == 0)
  
  for (t in seq.int(2L, numIter)) {
    
    if (t%%500 == 0){
      print(t)
    }
    
    ## consider each protein in turn
    for (i in seq.int(nrow(X))) {
      
      ## if assigned to a cluster remove statistics
      if ( Outlier[i, t - 1] == 1) {
        idx <- Component[i, t - 1] ## temporary variable index
        tempS <- mk[idx, ] %*% t(mk[idx, ]) ## temporary scatter, since mk to be update on next line
        mk[idx, ] <- (lambdak[idx] * mk[idx, ] - X[i, ]) / (lambdak[idx] - 1)
        lambdak[idx] <- lambdak[idx] - 1
        nuk[idx] <- nuk[idx] - 1
        nk[idx] <- nk[idx] - 1
        tau1 <- tau1 - 1
        sk[idx, , ] <- sk[idx, , ] - (X[i, ] %*% t(X[i, ])) +
          (lambdak[idx] + 1) * tempS  - lambdak[idx] * mk[idx,] %*% t(mk[idx,])
      } else {
        if (t > 2) {
          tau2 <- tau2 - 1
        }
      }
      
      ## compute probability of belonging to each organelle
      ## precompute terms for speed
      weight <- (nk + betak)/(sum(nk) + sum(betak) - 1) ## Component weights
      sigmak <- ((1 + lambdak) * sk)/(lambdak * (nuk - D + 1)) ## scale matrix
      degf <- nuk - D + 1 ## degrees freedom
      for (j in seq.int(K)) {
        ComponentProb[i, t, j] <- log(weight[j]) + pRoloc:::dmvtCpp(X[i, ,drop = FALSE],
                                                                    mu_ = mk[j, ],
                                                                    sigma_ = sigmak[j, , ],
                                                                    df_ = degf[j],
                                                                    log_ = TRUE,
                                                                    ncores_ = 1,
                                                                    isChol_ = FALSE)
      }
      
      ## normalise with underflow correction
      c <-  max(ComponentProb[i, t, ])
      ComponentProb[i, t , ] <- exp(ComponentProb[i ,t , ] - c) / sum(exp(ComponentProb[i, t, ] - c))
      
      ## sample component
      Component[i, t] <- sample(x = 1:K, size = 1, prob = ComponentProb[i, t , ] )
      
      ## compute outlier allocation
      n <- nrow(object)
      idk <- Component[i, t] ## temporary allocation variable
      OutlierProb[i, t, 1] <- log((tau1 + v)/(n + u + v - 1)) + pRoloc:::dmvtCpp(X[i, ,drop=FALSE],
                                                                                 mu_ = mk[idk, ],
                                                                                 sigma_ = sigmak[idk,,],
                                                                                 df_ = degf[idk],
                                                                                 log_ = TRUE,
                                                                                 ncores_ = 1,
                                                                                 isChol_ = FALSE)
      OutlierProb[i, t, 2] <- log((tau2 + u)/(n + u + v - 1)) + pRoloc:::dmvtCpp(X[i, ,drop = FALSE],
                                                                                 mu_ = M,
                                                                                 sigma_ = V,
                                                                                 df_ = 4,
                                                                                 log_ = TRUE,
                                                                                 ncores_ = 1,
                                                                                 isChol_ = FALSE)
      ## normalise and sample
      OutlierProb[i, t, ] <- exp(OutlierProb[i, t, ])/sum(exp(OutlierProb[i, t, ]))
      Outlier[i, t] <- sample(x = c(1, 0), 1, prob = OutlierProb[i, t, ]) ## reversed sample so 2nd entry is prob of 0.
      
      ## allocation statistics
      if ( Outlier[i, t] == 1) {
        idx <- Component[i, t] ## temporary variable index
        tempS <- mk[idx, ] %*% t(mk[idx, ]) ## temporary scatter, since mk to be update on next line
        mk[idx, ] <- (lambdak[idx] * mk[idx, ] + X[i, ]) / (lambdak[idx] + 1)
        lambdak[idx] <- lambdak[idx] + 1
        nuk[idx] <- nuk[idx] + 1
        nk[idx] <- nk[idx] + 1
        tau1 <- tau1 + 1
        #if ( t == 2) {
        #  tau2 <- tau2 - 1 ## default on first round is phi = 0
        #}
        sk[idx, , ] <- sk[idx,,] + (X[i,] %*% t(X[i,])) +
          (lambdak[idx] - 1) * tempS  - lambdak[idx] * mk[idx,] %*% t(mk[idx,])
      } else {
        #if (t > 2) {
        tau2 <- tau2 + 1
        #}
      }
    } ## end loop over proteins
  } ## end iterations
  
  ## create names for objects
  rownames(mk) <- names(lambdak) <- names(nuk) <- dimnames(sk)[[1]] <- c(markers, paste0("Phenotype ", seq.int(K_bar)))
  colnames(mk) <- dimnames(sk)[[2]] <- dimnames(sk)[[3]] <- sampleNames(object)
  
  ## name storage
  p_names <- rownames(X)
  rownames(Component) <- rownames(ComponentProb) <- rownames(Outlier) <- rownames(OutlierProb) <- p_names
  dimnames(ComponentProb)[[3]] <- c(markers, paste0("Phenotype ", seq.int(K_bar)))
  
  ## save Component parameters
  .ComponentParam <- pRoloc:::.ComponentParam(K = K, D = D,
                                              mk = mk,
                                              lambdak = lambdak,
                                              nuk = nuk,
                                              sk = sk)
  ## apply thinning and burn-in
  .Component <- Component[, retained]
  .ComponentProb <- ComponentProb[, retained, ]
  .Outlier <- Outlier[, retained]
  .OutlierProb <- OutlierProb[, retained, ]
  
  ## make MCMCChain object
  .MCMCChain <- pRoloc:::.MCMCChain(n = length(retained),
                                    K = K,
                                    N = nrow(X),
                                    Component = .Component,
                                    ComponentProb = .ComponentProb,
                                    Outlier = .Outlier,
                                    OutlierProb = .OutlierProb,
                                    ComponentParam = .ComponentParam)
  
  return(.MCMCChain)
}
