##' These functions implement the T augmented Gaussian mixture (TAGM)
##' model for mass spectrometry-based spatial proteomics datasets
##' using Markov-chain Monte-Carlo (MCMC) for inference.
##'
##' The `tagmMcmcTrain` function generates the samples from the
##' posterior distributions (object or class `MCMCParams`) based on an
##' annotated quantitative spatial proteomics dataset (object of class
##' [`MSnbase::MSnSet`]). Both are then passed to the `tagmPredict`
##' function to predict the sub-cellular localisation of protein of
##' unknown localisation. See the *pRoloc-bayesian* vignette for
##' details and examples. In this implementation, if numerical instability
##' is detected in the covariance matrix of the data a small multiple of
##' the identity is added. A message is printed if this conditioning step
##' is performed.
##'
##' @title Localisation of proteins using the TAGM MCMC method
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
##' @param numChains The number of parrallel chains to be run. Default
##'     it 4.
##' @param BPPARAM Support for parallel processing using the
##'     `BiocParallel` infrastructure. When missing (default), the
##'     default registered `BiocParallelParam` parameters are
##'     used. Alternatively, one can pass a valid `BiocParallelParam`
##'     parameter instance: `SnowParam`, `MulticoreParam`,
##'     `DoparParam`, ... see the `BiocParallel` package for
##'     details.
##' @return `tagmMcmcTrain` returns an instance of class
##'     `MCMCParams`.
##' @md
##' @references *A Bayesian Mixture Modelling Approach For Spatial
##'     Proteomics* Oliver M Crook, Claire M Mulvey, Paul D. W. Kirk,
##'     Kathryn S Lilley, Laurent Gatto bioRxiv 282269; doi:
##'     https://doi.org/10.1101/282269
##' @seealso The [plotEllipse()] function can be used to visualise
##'     TAGM models on PCA plots with ellipses.
##' @rdname tagm-mcmc
tagmMcmcTrain <- function(object,
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
                          BPPARAM = BiocParallel::bpparam()) {

    ## get expression marker data
    markersubset <- markerMSnSet(object, fcol = fcol)
    markers <- getMarkerClasses(markersubset, fcol = fcol)
    mydata <- exprs(markersubset)
    X <- exprs(unknownMSnSet(object, fcol = fcol))

    ## get data dize
    N <- nrow(mydata)
    D <- ncol(mydata)
    K <- length(markers)

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
                     FUN = tagmMcmcChain,
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
                     BPPARAM = BPPARAM)

    ## Construct class MCMCChains
    .ans <- .MCMCChains(chains = .res)

    ## Construct class MCMCParams
    out <- .MCMCParams(method = "TAGM.MCMC",
                       chains = .ans,
                       priors = .priors,
                       summary = .MCMCSummary())

    return(out)
}


##' @param params An instance of class `MCMCParams`, as generated by
##'     [tagmMcmcTrain()].
##' @param probJoint A `logical(1)` indicating whether to return the
##'     joint probability matrix, i.e. the probability for all classes
##'     as a new `tagm.mcmc.joint` feature variable.
##' @param probOutlier A `logical(1)` indicating whether to return the
##'     probability of being an outlier as a new `tagm.mcmc.outlier`
##'     feature variable. A high value indicates that the protein is
##'     unlikely to belong to any annotated class (and is hence
##'     considered an outlier).
##' @return `tagmMcmcPredict` returns an instance of class
##'     [`MSnbase::MSnSet`] containing the localisation predictions as
##'     a new `tagm.mcmc.allocation` feature variable. The allocation
##'     probability is encoded as `tagm.mcmc.probability`
##'     (corresponding to the mean of the distribution
##'     probability). In additionm the upper and lower quantiles of
##'     the allocation probability distribution are available as
##'     `tagm.mcmc.probability.lowerquantile` and
##'     `tagm.mcmc.probability.upperquantile` feature variables. The
##'     Shannon entropy is available in the `tagm.mcmc.mean.shannon`
##'     feature variable, measuring the uncertainty in the allocations
##'     (a high value representing high uncertainty; the highest value
##'     is the natural logarithm of the number of classes).
##' @md
##' @rdname tagm-mcmc
tagmMcmcPredict <- function(object,
                            params,
                            fcol = "markers",
                            probJoint = FALSE,
                            probOutlier = TRUE) {
    stopifnot(inherits(params, "MCMCParams"))
    ## Checks for object and MCMCParams match
    stopifnot(featureNames(unknownMSnSet(object, fcol = fcol))
              == rownames(params@summary@posteriorEstimates))

    ## Create marker set and size
    markerSet <- markerMSnSet(object, fcol = fcol)
    markers <- getMarkerClasses(object, fcol = fcol)
    M <- nrow(markerSet)
    K <- chains(params)[[1]]@K


    ## Get Summary object from MCMCParams maybe better to check
    ## columns exist/pass which objects we need
    .tagm.allocation <- c(as.character(params@summary@posteriorEstimates[,"tagm.allocation"]),
                          as.character(fData(markerSet)[, fcol]))
    .tagm.probability <- c(params@summary@posteriorEstimates[,"tagm.probability"],
                           rep(1, M)) ## set all probabilities of markers to 1.
    .tagm.probability.lowerquantile <- c(params@summary@posteriorEstimates[,"tagm.probability.lowerquantile"],
                                         rep(1, M)) ## set all probabilities of markers to 1.
    .tagm.probability.upperquantile <- c(params@summary@posteriorEstimates[,"tagm.probability.upperquantile"],
                                         rep(1, M)) ## set all probabilities of markers to 1.
    .tagm.mean.shannon <- c(params@summary@posteriorEstimates[,"tagm.mean.shannon"],
                            rep(0, M)) ## set all probabilities of markers to 0

    ## Create data frame to store new summaries
    .tagm.summary <- data.frame(tagm.mcmc.allocation = .tagm.allocation ,
                                tagm.mcmc.probability = .tagm.probability,
                                tagm.mcmc.probability.lowerquantile = .tagm.probability.lowerquantile,
                                tagm.mcmc.probability.upperquantile = .tagm.probability.upperquantile,
                                tagm.mcmc.mean.shannon = .tagm.mean.shannon)
    if (probOutlier)
        .tagm.summary$tagm.mcmc.outlier <-
            c(params@summary@posteriorEstimates[, "tagm.probability.Outlier"],
              rep(0, M)) ## set all probabilities of markers to 0


    ## Check number of rows match and add feature names
    stopifnot(nrow(.tagm.summary) == nrow(object))
    rownames(.tagm.summary) <- c(rownames(params@summary@posteriorEstimates),
                                 rownames(markerSet))

    ## Append data to fData of MSnSet
    fData(object) <- cbind(fData(object), .tagm.summary[rownames(fData(object)),])

    if  (probJoint) {
        ## create allocation matrix for markers
        .probmat <- matrix(0, nrow = nrow(markerSet), ncol = K)
        .class <- fData(markerSet)[, fcol]
        for (j in seq_len(nrow(markerSet))) {
            ## give markers prob 1
            .probmat[j, as.numeric(factor(.class), seq(1, length(unique(.class))))[j]] <- 1
        }
        colnames(.probmat) <- markers
        rownames(.probmat) <- rownames(markerSet)
        .joint <- rbind(params@summary@tagm.joint, .probmat)
        fData(object)$tagm.mcmc.joint <- .joint[rownames(fData(object)), ]
    }

    return(object)
}


##' @rdname tagm-mcmc
tagmPredict <- function(object,
                        params,
                        fcol = "markers",
                        probJoint = FALSE,
                        probOutlier = TRUE) {
    if (inherits(params, "MAPParams")) {
        ans <- tagmMapPredict(object, params, fcol, probJoint, probOutlier)
        ans@processingData@processing <- c(processingData(ans)@processing,
                                           paste0("Performed TAGM-MAP prediction", date()))
    }
    else if (inherits(params, "MCMCParams")) {
        ans  <- tagmMcmcPredict(object, params, fcol, probJoint, probOutlier)
        ans@processingData@processing <- c(processingData(ans)@processing,
                                           paste0("Performed TAGM-MCMC prediction", date()))
    } else
        stop("Parameters must either be 'MAPParams' or 'MCMCParams'.")
    return(ans)
}

##' @return `tagmMcmcProcess` returns an instance of class
##'     `MCMCParams` with its summary slot populated.
##' @md
##' @rdname tagm-mcmc
tagmMcmcProcess <- function(params) {
    ## get require slots
    myChain <- chains(params)[[1]]
    numChains <- length(chains(params))
    N <- myChain@N
    K <- myChain@K
    n <- myChain@n

    ## storage
    meanComponentProb = vector("list", length = numChains)
    meanOutlierProb = vector("list", length = numChains)
    tagm.joint <- matrix(0, nrow = N, ncol = K)
    tagm.outlier <- matrix(0, nrow = N, ncol = 2)
    .tagm.quantiles <- array(0, c(N, K, 2))
    pooled.Component <- matrix(0, nrow = N, ncol = n * numChains)
    pooled.ComponentProb <- array(0, c(N, n * numChains, K ))
    pooled.Outlier <- matrix(0, nrow = N, ncol = n * numChains)
    pooled.OutlierProb <- array(0, c(N, n * numChains, 2 ))

                                        # Calculate basic quantities
    for (j in seq_len(numChains)) {

        mc <- chains(params)[[j]]
        meanComponentProb[[j]] <- t(apply(mc@ComponentProb[, , ], 1, colMeans))
        meanOutlierProb[[j]] <- t(apply(mc@OutlierProb[, , ], 1, colMeans))
        tagm.joint <- tagm.joint + meanComponentProb[[j]]
        tagm.outlier <- tagm.outlier + meanOutlierProb[[j]]
        ## Pool chains
        pooled.Component[, n * (j - 1) + seq.int(n)] <- mc@Component
        pooled.ComponentProb[, n * (j - 1) + seq.int(n), ] <- mc@ComponentProb
        pooled.Outlier[, n * (j - 1)+ seq.int(n)] <- mc@Outlier
        pooled.OutlierProb[, n * (j - 1) + seq.int(n), ] <- mc@OutlierProb

    }
    ## take means across chains
    tagm.joint <- tagm.joint/numChains
    tagm.outlier <- tagm.outlier/numChains
    tagm.probability <- apply(tagm.joint, 1, max)
    tagm.allocation <- colnames(tagm.joint)[apply(tagm.joint, 1, which.max)]

    ## Calculate quantiles
    for (i in seq_len(N)) {
        for (j in seq_len(K)) {
            .tagm.quantiles[i, j, ] <- quantile(pooled.ComponentProb[i, , j], probs = c(0.025, 0.975))
        }
    }

    ## Store quantiles
    tagm.probability.lowerquantile <- .tagm.quantiles[cbind(1:N, apply(tagm.joint, 1, which.max), rep(1, N))]
    tagm.probability.upperquantile <- .tagm.quantiles[cbind(1:N, apply(tagm.joint, 1, which.max), rep(2, N))]

    ## Compute Shannon Entropy
    tagm.shannon <- -apply(pooled.ComponentProb * log(pooled.ComponentProb), c(1,2), sum)
    tagm.shannon[is.na(tagm.shannon)] <- 0
    tagm.mean.shannon <- rowMeans(tagm.shannon)

    ## Name entries
    names(tagm.mean.shannon) <- names(tagm.probability.lowerquantile) <-
        names(tagm.probability.upperquantile) <- rownames(myChain@Component)
    rownames(tagm.joint) <- names(tagm.probability) <-
        names(tagm.allocation) <- rownames(myChain@Component)
    colnames(tagm.joint) <- rownames(myChain@ComponentParam@mk)

    ## Outlier colNames
    colnames(tagm.outlier) <- c("tagm.probability.notOutlier", "tagm.probability.Outlier")

    ## Compute convergence diagnostics
    outlierTotal <- vector("list", length = numChains)
    for(j in seq_len(numChains)) {
        mc <- chains(params)[[j]]
        outlierTotal[[j]] <- coda::mcmc(colSums(mc@Outlier))
    }

    ## Summary of posterior estimates
    tagm.summary <- data.frame(tagm.allocation, tagm.probability, tagm.outlier,
                               tagm.probability.lowerquantile,
                               tagm.probability.upperquantile, tagm.mean.shannon)

    ## Compute covergence diagnostics
    .diagnostics <- matrix(0, nrow = 1, ncol = 2)
    if(numChains > 1){
        outlierTotal <- coda::as.mcmc.list(outlierTotal)
        gd <- coda::gelman.diag(x = outlierTotal, autoburnin = F)
        .diagnostics <- gd$psrf
        rownames(.diagnostics) <- c("outlierTotal")
    }

    ## Constructor for summary
    params@summary <- .MCMCSummary(posteriorEstimates = tagm.summary,
                                   diagnostics = .diagnostics,
                                   tagm.joint = tagm.joint)

    return(params)
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
##' @return `tagmMcmcChain` returns an instance of class
##'     [MCMCChain()].
##' @md
##' @noRd
tagmMcmcChain <- function(object,
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
                          v = 10) {
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
    K <- length(markers)

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
    nk <- c(table(fData(markersubset)[, fcol])[markers])

    for (j in seq.int(K))
        xk[j, ] <- colSums(mydata[fData(markersubset)[, fcol] == markers[j], ])/nk[j]

    lambdak <- lambda0 + nk
    nuk <- nu0 + nk
    mk <- t((t(nk * xk) + lambda0 * mu0)) / lambdak
    betak <- beta0 + nk

    for(j in seq.int(K))
        sk[j, , ] <- S0 + t(mydata[fData(markersubset)[, fcol] == markers[j], ]) %*%
            mydata[fData(markersubset)[, fcol] == markers[j],] +
            lambda0 * mu0 %*% t(mu0) - lambdak[j] * mk[j, ] %*% t(mk[j, ])


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
        ComponentProb[, 1, j] <- dmvtCpp(X,
                                         mu_ = mk[j, ],
                                         sigma_ = (1 + lambdak[j]) * sk[j, , ] / (lambdak[j] * (nuk[j] - D + 1)),
                                         df_ = nuk[j] - D + 1,
                                         log_ = TRUE,
                                         ncores_ = 1,
                                         isChol_ = FALSE)

    Component[, 1] <- apply(X = ComponentProb[, 1, ], 1, FUN = which.max)

    ## initially assign all components to global component i.e. allocstruc = 0
    Outlier[, 1] <- 0

    ## initial allocation statistics
    tau1 <- sum(Outlier[, 1] == 1) + N
    tau2 <- sum(Outlier[, 1] == 0)

    for (t in seq.int(2L, numIter)) {

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
                ComponentProb[i, t, j] <- log(weight[j]) + dmvtCpp(X[i, ,drop = FALSE],
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
            OutlierProb[i, t, 1] <- log((tau1 + v)/(n + u + v - 1)) + dmvtCpp(X[i, ,drop=FALSE],
                                                                              mu_ = mk[idk, ],
                                                                              sigma_ = sigmak[idk,,],
                                                                              df_ = degf[idk],
                                                                              log_ = TRUE,
                                                                              ncores_ = 1,
                                                                              isChol_ = FALSE)
            OutlierProb[i, t, 2] <- log((tau2 + u)/(n + u + v - 1)) + dmvtCpp(X[i, ,drop = FALSE],
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
                if ( t == 2) {
                    tau2 <- tau2 - 1 ## default on first round is phi = 0
                }
                sk[idx, , ] <- sk[idx,,] + (X[i,] %*% t(X[i,])) +
                    (lambdak[idx] - 1) * tempS  - lambdak[idx] * mk[idx,] %*% t(mk[idx,])
            } else {
                if (t > 2) {
                    tau2 <- tau2 + 1
                }
            }
        } ## end loop over proteins
    } ## end iterations

    ## create names for objects
    rownames(mk) <- names(lambdak) <- names(nuk) <- dimnames(sk)[[1]] <- markers
    colnames(mk) <- dimnames(sk)[[2]] <- dimnames(sk)[[3]] <- sampleNames(object)

    ## name storage
    p_names <- rownames(X)
    rownames(Component) <- rownames(ComponentProb) <- rownames(Outlier) <- rownames(OutlierProb) <- p_names
    dimnames(ComponentProb)[[3]] <- markers

    ## save Component parameters
    .ComponentParam <- .ComponentParam(K = K, D = D,
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
    .MCMCChain <- .MCMCChain(n = length(retained),
                             K = K,
                             N = nrow(X),
                             Component = .Component,
                             ComponentProb = .ComponentProb,
                             Outlier = .Outlier,
                             OutlierProb = .OutlierProb,
                             ComponentParam = .ComponentParam)

    return(.MCMCChain)
}
