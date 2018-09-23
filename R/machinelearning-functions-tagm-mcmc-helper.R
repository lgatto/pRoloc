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

mcmc_plot_probs <- function(x, fname, n = 1) {
    stopifnot(inherits(x, "MCMCParams"))
    stopifnot(require("ggplot2"))
    chain <- pRoloc:::chains(x)[[n]]
    dfr <- as.data.frame(chain@ComponentProb[fname, , ])
    colnames(dfr) <- rownames(chain@ComponentParam@mk)
    dfr_long <- data.frame(Organelle = rep(names(dfr), each = nrow(dfr)),
                           Probability = unlist(dfr, use.names = FALSE),
                           row.names = NULL,
                           stringsAsFactors = FALSE)
    gg2 <- ggplot(dfr_long,
                  aes(Organelle, Probability,
                      width = (Probability))) +
        geom_violin(aes(fill = Organelle), scale = "width")
    gg2 <- gg2 +
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
