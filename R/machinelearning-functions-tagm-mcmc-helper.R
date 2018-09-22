mcmc_get_outliers <- function(x) {
    stopifnot(inherits(x, "MCMCParams"))
    lapply(x@chains@chains, function(mc) coda::mcmc(colSums(mc@Outlier)))
}

mcmc_thin_chains <- function(x, n) {
    stopifnot(inherits(x, "MCMCParams"))
    .chain <- pRoloc:::chains(newTanMcmc)[[1]]
    K <- .chain@K # Number of components
    N <- .chain@N # Number of Proteins
    chainlist <-
        lapply(x@chains@chains, function(chain) {
            .ComponentParam <- chain@ComponentParam
            ## Subset MCMC iterations
            retain <- seq.int(n + 1, chain@n) ## samples to retain
            ## Check correct number of iterations
            stopifnot(ncol(chain@Component[, retain]) == n)
            .Component <- chain@Component[, retain]
            .ComponentProb <- chain@ComponentProb[, retain, ]
            .Outlier <- chain@Outlier[, retain]
            .OutlierProb <- chain@OutlierProb[, retain, ]
            pRoloc:::.MCMCChain(n = as.integer(n),
                                K = K,
                                N = N,
                                Component = .Component,
                                ComponentProb = .ComponentProb,
                                Outlier = .Outlier,
                                OutlierProb = .OutlierProb,
                                ComponentParam = .ComponentParam)
        })

    mcmc_chainlist <- pRoloc:::.MCMCChains(chains = chainlist)
    pRoloc:::.MCMCParams(method = "TAGM.MCMC",
                         chains = mcmc_chainlist,
                         priors = x@priors,
                         summary = pRoloc:::.MCMCSummary())
}
