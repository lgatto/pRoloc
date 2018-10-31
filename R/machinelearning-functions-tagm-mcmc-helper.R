
# Helper function to get number of outlier at each iteration of MCMC
mcmc_get_outliers <- function(x) {
    stopifnot(inherits(x, "MCMCParams"))
    lapply(x@chains@chains, function(mc) coda::mcmc(colSums(1 - mc@Outlier)))
}

# Helper function to get mean component allocation at each iteration of MCMC

mcmc_get_meanComponent <- function(x) {
  stopifnot(inherits(x, "MCMCParams"))
  lapply(x@chains@chains, function(mc) coda::mcmc(colMeans(mc@Component)))
}

# Helper function to get mean probability of belonging to outlier at each iteration

mcmc_get_meanoutliersProb <- function(x) {
  stopifnot(inherits(x, "MCMCParams"))
  lapply(x@chains@chains, function(mc) coda::mcmc(colMeans(mc@OutlierProb[, ,2])))
}

# Wrapper for the geweke diagnostics from coda package also return p-values

geweke_test <- function(x) {
  res <- matrix(NA, nrow = 2, ncol = length(x))
  gwk <- sapply(x, coda::geweke.diag, simplify = TRUE)
  res[1, ] <- unlist(gwk[1, ])
  res[2, ] <- pnorm(abs(unlist(gwk[1, ])), lower.tail=FALSE) * 2
  colnames(res) <- paste0("chain ", seq.int(x))
  rownames(res) <- c("z.value", "p.value")
  return(res)
}

# Helper function to pool chains together after processing

mcmc_pool_chains <- function(param) {
  stopifnot(inherits(param, "MCMCParams"))
  
  param1 <- pRoloc:::chains(param)[[1]]

  n <- param1@n
  nPool <- param1@n * length(param) # total number of iteration increase
  KPool <- param1@K  # number of components unchanged
  NPool <- param1@N # number of proteins doesn't change
  numChains <- length(param)
  
  pooled.Component <- matrix(0, nrow = NPool, ncol = nPool)
  pooled.ComponentProb <- array(0, c(NPool, nPool, KPool ))
  pooled.Outlier <- matrix(0, nrow = NPool, ncol = nPool)
  pooled.OutlierProb <- array(0, c(NPool, nPool, 2 ))
  
  rownames(pooled.Component) <- rownames(param1@Component)
  rownames(pooled.ComponentProb) <- rownames(param1@ComponentProb)
  rownames(pooled.Outlier) <- rownames(param1@Outlier)
  rownames(pooled.OutlierProb) <- rownames(param1@OutlierProb)
  dimnames(pooled.ComponentProb)[[3]] <- dimnames(param1@ComponentProb)[[3]]
  
  
  # Calculate basic quantities
  for (j in seq_len(numChains)) {
    
    mc <- pRoloc:::chains(param)[[j]]
    ## Pool chains
    pooled.Component[, n * (j - 1) + seq.int(n)] <- mc@Component
    pooled.ComponentProb[, n * (j - 1) + seq.int(n), ] <- mc@ComponentProb
    pooled.Outlier[, n * (j - 1)+ seq.int(n)] <- mc@Outlier
    pooled.OutlierProb[, n * (j - 1) + seq.int(n), ] <- mc@OutlierProb
    
  }
  
  mk.list <- lapply(param@chains@chains,function(x) x@ComponentParam@mk)
  lambdak.list <- lapply(param@chains@chains,function(x) x@ComponentParam@lambdak)
  nuk.list <- lapply(param@chains@chains,function(x) x@ComponentParam@nuk)
  sk.list <- lapply(param@chains@chains,function(x) x@ComponentParam@sk)
  
  ## save Component parameters
  .ComponentParam <- pRoloc:::.ComponentParam(K = KPool, D = param1@ComponentParam@D,
                                     mk = Reduce("+", mk.list) / length(mk.list),
                                     lambdak = Reduce("+", lambdak.list) / length(lambdak.list),
                                     nuk = Reduce("+", nuk.list) / length(nuk.list),
                                     sk = Reduce("+", sk.list) / length(sk.list))
  ## apply thinning and burn-in
  .Component <- pooled.Component
  .ComponentProb <- pooled.ComponentProb
  .Outlier <- pooled.Outlier
  .OutlierProb <- pooled.OutlierProb
  
  ## make MCMCChain object
  .MCMCChain <- pRoloc:::.MCMCChain(n = nPool,
                           K = KPool,
                           N = NPool,
                           Component = .Component,
                           ComponentProb = .ComponentProb,
                           Outlier = .Outlier,
                           OutlierProb = .OutlierProb,
                           ComponentParam = .ComponentParam)
  
  ## Make MCMCChains with single object
  .MCMCChains <- pRoloc:::.MCMCChains(chains = list(.MCMCChain))
  
  ## Make MCMCParams object
  pRoloc:::.MCMCParams(method = "TAGM.MCMC",
                       chains = .MCMCChains,
                       priors = param@priors,
                       summary = pRoloc:::.MCMCSummary())
  
}

# Helper function to burn n iterations from the front of the chains

mcmc_burn_chains <- function(x, n) {
    stopifnot(inherits(x, "MCMCParams"))
    .chain <- pRoloc:::chains(x)[[1]]
    K <- .chain@K # Number of components
    N <- .chain@N # Number of Proteins
    chainlist <-
        lapply(x@chains@chains, function(chain) {
            .ComponentParam <- chain@ComponentParam
            ## Subset MCMC iterations
            retain <- seq.int(n + 1, chain@n) ## samples to retain
            ## Check correct number of iterations
            stopifnot(ncol(chain@Component[, retain]) == (chain@n - n))
            .Component <- chain@Component[, retain]
            .ComponentProb <- chain@ComponentProb[, retain, ]
            .Outlier <- chain@Outlier[, retain]
            .OutlierProb <- chain@OutlierProb[, retain, ]
            pRoloc:::.MCMCChain(n = as.integer(chain@n - n),
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

# Helper function to subsample the chains, known informally as thinning.

mcmc_thin_chains <- function(x, freq) {
  stopifnot(inherits(x, "MCMCParams"))
  .chain <- pRoloc:::chains(x)[[1]]
  K <- .chain@K # Number of components
  N <- .chain@N # Number of Proteins
  nThin <- floor(.chain@n/freq) # Number of iterations after thinning
  chainlist <-
    lapply(x@chains@chains, function(chain) {
      .ComponentParam <- chain@ComponentParam
      ## Subset MCMC iterations
      retain <- freq * seq.int(1:nThin) ## samples to retain
      ## Check correct number of iterations
      stopifnot(ncol(chain@Component[, retain]) == nThin)
      .Component <- chain@Component[, retain]
      .ComponentProb <- chain@ComponentProb[, retain, ]
      .Outlier <- chain@Outlier[, retain]
      .OutlierProb <- chain@OutlierProb[, retain, ]
      pRoloc:::.MCMCChain(n = as.integer(nThin),
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
                       priors = param@priors,
                       summary = pRoloc:::.MCMCSummary())
}

# Plotting method for violins

setMethod("plot", c("MCMCParams", "character"),
          function(x, y, ...) {
              mcmc_plot_probs(x, y, n = 1)
          })
# Plotting function for violins using ggplot2.

mcmc_plot_probs <- function(param, fname, n = 1) {
    stopifnot(inherits(param, "MCMCParams"))
    stopifnot(length(fname) == 1)
    chain <- pRoloc:::chains(param)[[n]]
    stopifnot(fname %in% rownames(chain@ComponentProb))
    dfr <- as.data.frame(chain@ComponentProb[fname, , ])
    colnames(dfr) <- rownames(chain@ComponentParam@mk)
    dfr_long <- data.frame(Organelle = rep(names(dfr), each = nrow(dfr)),
                           Probability = unlist(dfr, use.names = FALSE),
                           row.names = NULL,
                           stringsAsFactors = FALSE)
    gg2 <- ggplot(dfr_long,
                  aes(Organelle, Probability,
                      width = (Probability))) +
        geom_violin(aes(fill = Organelle), scale = "width", bw = 0.05)
    gg2 <- gg2 + theme_bw() + 
        scale_fill_manual(values = pRoloc::getStockcol()[seq_len(nrow(dfr))]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.title.x = element_blank()) 
    gg2 <- gg2 +
        ylab("Membership Probability") +
        ggtitle(paste0("Distribution of Subcellular Membership for Protein ", fname ))
    gg2 <- gg2 +
        theme(legend.position = "none")
    return(gg2)
}
