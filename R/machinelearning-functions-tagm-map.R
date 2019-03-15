##' These functions implement the T augmented Gaussian mixture (TAGM)
##' model for mass spectrometry-based spatial proteomics datasets
##' using the maximum a posteriori (MAP) optimisation routine.
##'
##' The `tagmMapTrain` function generates the MAP parameters (object or class
##' `MAPParams`) based on an annotated quantitative spatial proteomics dataset
##' (object of class [`MSnbase::MSnSet`]). Both are then passed to the
##' `tagmPredict` function to predict the sub-cellular localisation of protein
##' of unknown localisation. See the *pRoloc-bayesian* vignette for details and
##' examples. In this implementation, if numerical instability is detected in
##' the covariance matrix of the data a small multiple of the identity is
##' added. A message is printed if this conditioning step is performed.
##'
##' @title Localisation of proteins using the TAGM MAP method
##'
##' @param object An [`MSnbase::MSnSet`] containing the spatial
##'     proteomics data to be passed to `tagmMapTrain` and
##'     `tagmPredict`.
##' @param fcol The feature meta-data containing marker definitions.
##'     Default is `markers`.
##' @param method A `charachter()` describing the inference method for
##'     the TAGM algorithm. Default is `"MAP"`.
##' @param numIter The number of iterations of the
##'     expectation-maximisation algorithm. Default is 100.
##' @param mu0 The prior mean. Default is `colMeans` of the expression
##'     data.
##' @param lambda0 The prior shrinkage. Default is 0.01.
##' @param nu0 The prior degreed of freedom. Default is
##'     `ncol(exprs(object)) + 2`
##' @param S0 The prior inverse-wishary scale matrix. Empirical prior
##'     used by default.
##' @param beta0 The prior Dirichlet distribution
##'     concentration. Default is 1 for each class.
##' @param u The prior shape parameter for Beta(u, v). Default is 2
##' @param v The prior shape parameter for Beta(u, v). Default is 10.
##' @param seed The optional random number generator seed.
##' @return `tagmMapTrain` returns an instance of class [MAPParams()].
##' @md
##' @references *A Bayesian Mixture Modelling Approach For Spatial
##'     Proteomics* Oliver M Crook, Claire M Mulvey, Paul D. W. Kirk,
##'     Kathryn S Lilley, Laurent Gatto bioRxiv 282269; doi:
##'     https://doi.org/10.1101/282269
##' @author Oliver M. Crook
##' @seealso The [plotEllipse()] function can be used to visualise
##'     TAGM models on PCA plots with ellipses. The [tagmMapTrain()]
##'     function to use the TAGM MAP method.
##' @rdname tagm-map
tagmMapTrain <- function(object,
                         fcol = "markers",
                         method = "MAP",
                         numIter = 100,
                         mu0 = NULL,
                         lambda0 = 0.01,
                         nu0 = NULL,
                         S0 = NULL,
                         beta0 = NULL,
                         u = 2,
                         v = 10,
                         seed = NULL) {

    ## get expression marker data
    markersubset <- markerMSnSet(object, fcol = fcol)
    markers <- getMarkerClasses(markersubset, fcol = fcol)
    mydata <- exprs(markersubset)
    X <- exprs(unknownMSnSet(object, fcol = fcol))

    if (is.null(seed))
        seed <- sample(.Machine$integer.max, 1)
    .seed <- as.integer(seed)
    set.seed(.seed)


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
    muk <- matrix(0, nrow = K, ncol = D)
    sigmak <- array(0, dim = c(K, D, D))
    xk <- matrix(0, nrow = K, ncol = D)

    ## update prior with training data
    nk <- tabulate(fData(markersubset)[, fcol])

    for (j in seq.int(K))
        xk[j, ] <- colSums(mydata[fData(markersubset)[, fcol] == markers[j], ])/nk[j]

    lambdak <- lambda0 + nk
    nuk <- nu0 + nk
    mk <- t((t(nk * xk) + lambda0 * mu0)) / lambdak

    for(j in seq.int(K))
        sk[j, , ] <- S0 + t(mydata[fData(markersubset)[, fcol] == markers[j], ]) %*%
            mydata[fData(markersubset)[, fcol] == markers[j],] +
            lambda0 * mu0 %*% t(mu0) - lambdak[j] * mk[j, ] %*% t(mk[j, ])

    betak <- beta0 + nk

    ## initial posterior mode
    muk <- mk
    for (j in seq.int(K))
        sigmak[j, , ] <- sk[j, , ] / (nuk[j] + D + 1)

    ## initial cluster probabilty weights
    pik <- (betak - 1) / (sum(betak) - K)

    ## global parameters
    M <- colMeans(exprs(object))
    V <- cov(exprs(object))/2
    eigsV <- eigen(V)
    if (min(eigsV$values) < .Machine$double.eps) {
      V <- cov(exprs(object))/2 + diag(10^{-6}, D)
      message("co-linearity detected; a small multiple of
              the identity was added to the covariance")
    }
    eps <- (u - 1) / (u + v - 2)

    ## storage for Estep
    a <- matrix(0, nrow = nrow(X), ncol = K)
    b <- matrix(0, nrow = nrow(X), ncol = K)
    w <- matrix(0, nrow = nrow(X), ncol = K)
    ## storage for Mstep
    xbar <- matrix(0, nrow = nrow(X), ncol = K)
    lambda <- matrix(0, K)
    nu <- matrix(0, K)
    m <- matrix(0, K, D)
    S <- array(0, c(K, D, D))
    loglike <- vector(mode = "numeric", length = numIter)

    for (t in seq.int(numIter)) {
        ## E-Step, log computation to avoid underflow
        for (k in seq.int(K)) {
            a[, k] <- log( pik[k] ) +
                log( 1 - eps) +
                mvtnorm::dmvnorm(X,
                                 mean = muk[k, ],
                                 sigma = sigmak[k, , ],
                                 log = TRUE)
            b[, k] <- log( pik[k] ) +
                log(eps) +
                mvtnorm::dmvt(X,
                              delta = M,
                              sigma = V,
                              df = 4,
                              log = TRUE)
        }

        ## correct for underflow by adding constant
        ab <- cbind(a,b)
        c <- apply(ab, 1, max)
        ab <- ab - c                   # add constant
        ab <- exp(ab)/rowSums(exp(ab)) # normlise
        a <- ab[, 1:K]
        b <- ab[, (K + 1):(2 * K)]
        w <- a + b
        r <- colSums(w)

        ## M-Step
        ## structure weights
        eps <- (u + sum(b) - 1) / ( (sum(a) + sum(b)) + (u + v) - 2)
        xbar <- apply(a, 2, function(x) colSums( x * X ))
        xbar[, colSums(xbar) != 0] <- t(t(xbar[, colSums(xbar)!=0])/colSums(a)[colSums(xbar)!=0])

        ## component weights
        pik <- (r + betak - 1)/(nrow(X) + sum(betak) - K)

        ## component parameters
        lambda <- lambdak + colSums(a)
        nu <- nuk + colSums(a)
        m <- (colSums(a) * t(xbar) + lambdak * mk) / lambda

        ## comptute scatter matrix
        TS <- array(0, c(K, D, D)) # temporary storage
        for(j in seq.int(K)) {
            for(i in seq.int(nrow(X))) {
                TS[j, , ] <- TS[j, , ] + a[i, j] * (X[i, ] - xbar[, j]) %*% t((X[i, ] - xbar[, j]))
            }
        }

        ## compute variance-covariance parameters
        vv <- (lambdak * colSums(a))/ lambda # temporary shrinkage variable
        for (j in seq.int(K)) {
            S[j, , ] <- sk[j, , ] + vv[j] * (xbar[, j] - mk[j,]) %*%
                t((xbar[, j] - mk[j,])) + TS[j, , ]
            sigmak[j, , ] <- S[j, , ]/(nu[j] + D + 2)
        }
        muk <- m

        ## compute log-likelihood, using recursive addition method
        for (j in seq.int(K)) {
            loglike[t] <- loglike[t] +
                sum( a[, j] * mvtnorm::dmvnorm(X, mean = muk[j, ], sigma = sigmak[j, , ], log = TRUE)) +
                sum( w[,j] * log(pik[j]) ) +
                LaplacesDemon::dinvwishart(Sigma = sigmak[j, , ], nu = nu0, S = S0, log = TRUE) +
                mvtnorm::dmvnorm(muk[j, ], mean = mu0, sigma = sigmak[j, , ], log = TRUE)

        }
        loglike[t] <- loglike[t] + sum(a) * log(1 - eps) + sum(b) * log(eps) +
            sum(rowSums(b) *  mvtnorm::dmvt(X, delta = M, sigma = V, df = 4, log = TRUE)) +
            stats::dbeta(eps, shape1 = u, shape2 = v, log = TRUE) +
            LaplacesDemon::ddirichlet(pik, alpha = beta0/K, log = TRUE)


    }

    ## save MAP estimates and log posterior
    .posteriors <- list(mu = muk,
                        sigma = sigmak,
                        weights = pik,
                        epsilon = eps,
                        logposterior = loglike)

    new("MAPParams",
        method = method,
        seed = .seed,
        priors = .priors,
        posteriors = .posteriors,
        datasize = list(
            "data" = dim(object)))
}



##' @param params An instance of class [`MAPParams`], as generated by
##'     [tagmMapTrain()].
##' @param probJoint A `logical(1)` indicating whether to return the
##'     joint probability matrix, i.e. the probability for all classes
##'     as a new `tagm.map.joint` feature variable.
##' @param probOutlier A `logical(1)` indicating whether to return the
##'     probability of being an outlier as a new `tagm.map.outlier`
##'     feature variable. A high value indicates that the protein is
##'     unlikely to belong to any annotated class (and is hence
##'     considered an outlier).
##' @return `tagmPredict` returns an instance of class
##'     [`MSnbase::MSnSet`] containing the localisation predictions as
##'     a new `tagm.map.allocation` feature variable.
##' @md
##' @rdname tagm-map
tagmMapPredict <- function(object,
                           params,
                           fcol = "markers",
                           probJoint = FALSE,
                           probOutlier = TRUE) {
    stopifnot(inherits(params, "MAPParams"))
    ## get parameters from
    posteriors <- params@posteriors
    eps <- posteriors$epsilon
    mu <- posteriors$mu
    sigma <- posteriors$sigma
    weights <- posteriors$weights

    ## split unknowns and markers
    unknownsubset <- unknownMSnSet(object, fcol = fcol)
    markersubset <- markerMSnSet(object, fcol = fcol)
    markers <- getMarkerClasses(markersubset, fcol = fcol)

    ## get data to predict
    X <- exprs(unknownsubset)
    K <- length(markers)
    D <- ncol(X)

    a <- matrix(0, nrow = nrow(X), ncol = K)
    b <- matrix(0, nrow = nrow(X), ncol = K)
    predictProb <- matrix(0, nrow = nrow(X), ncol =  K)
    organelleAlloc <- data.frame(pred = rep(NA_character_, nrow(X)),
                                 prob = rep(NA_real_, nrow(X)))

    ## global parameters
    M <- colMeans(exprs(object))
    V <- cov(exprs(object))/2
    eigsV <- eigen(V)
    if (min(eigsV$values) < .Machine$double.eps) {
      V <- cov(exprs(object))/2 + diag(10^{-6}, D)
    }

    for (j in seq.int(K)) {
        a[, j] <- log( weights[j] ) +
            log( 1 - eps) +
            mvtnorm::dmvnorm(X,
                             mean = mu[j, ],
                             sigma = sigma[j, , ],
                             log = TRUE)
        b[, j] <- log( weights[j] ) +
            log(eps) +
            mvtnorm::dmvt(X,
                          delta = M,
                          sigma = V,
                          df = 4,
                          log = TRUE)
    }

    ## correct for underflow by adding constant
    ab <- cbind(a, b)
    c <- apply(ab, 1, max)
    ab <- ab - c                   # add constant
    ab <- exp(ab)/rowSums(exp(ab)) # normlise
    a <- ab[, 1:K]
    b <- ab[, (K + 1):(2 * K)]
    .predictProb <- a + b

    colnames(.predictProb) <- markers

    organelleAlloc[["pred"]] <- markers[apply(a, 1, which.max)]
    probAlloc <- apply(a, 1, which.max)

    for (i in seq.int(nrow(X))) {
        organelleAlloc[i, "prob"] <- a[i, probAlloc[i]]
    }
    rownames(b) <- rownames(a) <- rownames(unknownsubset)
    rownames(organelleAlloc) <- rownames(.predictProb) <- rownames(unknownsubset)

    ## predicted classes and probabilities (markers set to 1)
    .pred <- c(organelleAlloc[["pred"]], as.character(fData(markersubset)[, fcol]))
    .prob <- c(organelleAlloc[["prob"]], rep(1, nrow(markersubset)))

    ## probablity of being an outlier (markers set to 0)
    .outlier <- c(rowSums(b), rep(0, nrow(markersubset)))

    ## making sure rownames align
    names(.outlier) <-
        names(.pred) <-
        names(.prob) <-
        c(rownames(unknownsubset), rownames(markersubset))

    fData(object)$tagm.map.allocation <- .pred[rownames(fData(object))]
    fData(object)$tagm.map.probability <- .prob[rownames(fData(object))]

    if  (probJoint) {
        ## create allocation matrix for markers
        .probmat <- matrix(0, nrow = nrow(markerMSnSet(object)), ncol = K)
        .class <- fData(markerMSnSet(object))[, fcol]
        for (j in seq_len(nrow(markerMSnSet(object)))) {
            ## give markers prob 1
            .probmat[j, as.numeric(factor(.class), seq(1,length(unique(.class))))[j]] <- 1
        }
        colnames(.probmat) <- markers
        rownames(.probmat) <- rownames(markersubset)
        .joint <- rbind(.predictProb, .probmat)
        fData(object)$tagm.map.joint <- .joint[rownames(fData(object)), ]
    }

    if (probOutlier)
        fData(object)$tagm.map.outlier <- .outlier[rownames(fData(object))]

    if (validObject(object))
        return(object)
}
