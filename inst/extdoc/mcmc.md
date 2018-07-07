
```r
library("pRoloc")
```

```
## Error: package 'MSnbase' could not be loaded
```

```r
library("pRolocdata")
```

```
## Error: package 'MSnbase' could not be loaded
```

```r
require("reshape")
require("ggplot2")
require("coda")

setStockcol(paste0(getStockcol(), 90))
```

```
## Error in setStockcol(paste0(getStockcol(), 90)): could not find function "setStockcol"
```

# *TAGM-MCMC*, in details

This section explains how to manually manipulate the MCMC output of
the Bayesian model TAGM-MCMC applied to spatial proteomics data
[@Crook:2018]. First, we load the MCMC data (as produced by
`tagmMcmcTrain`) to be used and the packages required for analysis.

The `tanTagm.rda` file is nearly 400MB large and isn't direcly
distributed in a package but can be downloaded from the
[http://bit.ly/tagm-mcmc](http://bit.ly/tagm-mcmc) google drive to
reproduce the following analyses. Alternatively, to avoid manual
intervention, the following code chunk downloads the file from an
alternative server:


```r
destdir <- tempdir()
destfile <- file.path(destdir, "tanTagm.rda")
download.file("http://proteome.sysbiol.cam.ac.uk/lgatto/files/tanTagm.rda",
			  destfile)
load(destfile)
```


```r
## Load tanTagm data containing the MCMC analysis, which is assumed to
## be available in the working directory.
load("tanTagm.rda")
tanTagm
```

```
## Loading required package: pRoloc
```

```
## Loading required package: MSnbase
```

```
## Error: package or namespace load failed for 'MSnbase' in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]):
##  there is no package called 'BiocInstaller'
```

```
## Error in .requirePackage(package): unable to find required package 'pRoloc'
```

We now load the example data for which we performed this Bayesian
analysis. This Drosophila embryo spatial proteomics experiment is
from [@Tan2009].


```r
data(tan2009r1) ## get data from pRolocdata
```

```
## Warning in data(tan2009r1): data set 'tan2009r1' not found
```

The `tanTagm` data was produce by executing the `tagmMcmcTrain`
function. $20000$ iterations were performed, automatically discarding
$5000$ iterations for burnin, sub-sampling the chains by $10$ and
running a total of $4$ chains in parallel. This results in $4$ chains
each with $1500$ MCMC iterations.


```r
tanTagm <- tagmMcmcTrain(object = tan2009r1,
						 numIter = 20000,
						 burnin = 5000,
						 thin = 10,
						 numChains = 4)
```

## Data exploration and convergence diagnostics

Using parallel chains for analysis allows us to diagnose whether our
MCMC algorithm has converged (or not). The number of parallel MCMC
chains used in this analysis was 4.


```r
## Get number of chains
nChains <- length(tanTagm@chains)
```

```
## Loading required package: pRoloc
```

```
## Loading required package: MSnbase
```

```
## Error: package or namespace load failed for 'MSnbase' in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]):
##  there is no package called 'BiocInstaller'
```

```
## Error in .requirePackage(package): unable to find required package 'pRoloc'
```

```r
nChains
```

```
## Error in eval(expr, envir, enclos): object 'nChains' not found
```

The following code chunks sets up a manual convegence diagnostic
check. We make use of objects and methods in the package
*[coda](https://CRAN.R-project.org/package=coda)* to peform this analysis [@coda].  We
calculate the total number of outliers at each iteration of each chain
and if the algorithm has converged this number should be the same (or
very similar) across all 4 chains. We can observe this by sight by
producing trace plots for each MCMC chain.


```r
## Convergence diagnostic to see if more we need to discard any
## iterations or entire chains.
outlierTotal <- vector("list", length = nChains)
```

```
## Error in vector("list", length = nChains): object 'nChains' not found
```

```r
## Compute the number of outliers for each iteration for each chain
for (j in seq_len(nChains)) {
  mc <- pRoloc:::chains(tanTagm)[[j]]
  outlierTotal[[j]] <- coda::mcmc(colSums(mc@Outlier))
}
```

```
## Error in eval(expr, envir, enclos): object 'nChains' not found
```

```r
## Carefully using coda S3 objects to produce trace plots and histograms
plot(outlierTotal[[1]], col = "blue", main = "Chain 1")
```

```
## Error in plot(outlierTotal[[1]], col = "blue", main = "Chain 1"): object 'outlierTotal' not found
```

```r
plot(outlierTotal[[2]], col = "red", main = "Chain 2")
```

```
## Error in plot(outlierTotal[[2]], col = "red", main = "Chain 2"): object 'outlierTotal' not found
```

```r
plot(outlierTotal[[3]], col = "green", main = "Chain 3")
```

```
## Error in plot(outlierTotal[[3]], col = "green", main = "Chain 3"): object 'outlierTotal' not found
```

```r
plot(outlierTotal[[4]], col = "orange", main = "Chain 4")
```

```
## Error in plot(outlierTotal[[4]], col = "orange", main = "Chain 4"): object 'outlierTotal' not found
```

We can use the *[coda](https://CRAN.R-project.org/package=coda)* package to produce
summaries of our chains. Here is the `coda` summary for the first
chain.


```r
## all chains average around 360 outliers
summary(outlierTotal[[1]])
```

```
## Error in summary(outlierTotal[[1]]): object 'outlierTotal' not found
```

In this case our chains looks very good. They all oscillate around an
average of 360 outliers. There is no observed monotonicity in our
output. However, for a more rigorous and unbiased analysis of
convergence we can calculate the Gelman diagnostics using the
*[coda](https://CRAN.R-project.org/package=coda)* package
[@Gelman:1992,@Brools:1998]. This statistics is often refered to as
$\hat{R}$ or the potential scale reduction factor. The idea of the
Gelman diagnostics is to so compare the inter and intra chain
variance. The ratio of these quantities should be close to one.  The
actual statistics computed is more complicated, but we do not go
deeper here and a more detailed and in depth discussion can be found
in the references. The *[coda](https://CRAN.R-project.org/package=coda)* package also
reports the $95\%$ upper confidence interval of the $\hat{R}$
statistic.


```r
## Can check gelman diagnostic for convergence (values less than <1.05
## are good for convergence)
gelman.diag(outlierTotal) # the Upper C.I. is 1 so mcmc has clearly converged
```

```
## Error in as.mcmc.list(x): object 'outlierTotal' not found
```

We can also look at the Gelman diagnostics statistics for pairs of chains.


```r
## We can also check individual pairs of chains for convergence
gelman.diag(outlierTotal[1:2]) # the upper C.I is 1.01
```

```
## Error in as.mcmc.list(x): object 'outlierTotal' not found
```

```r
gelman.diag(outlierTotal[c(1,3)]) # the upper C.I is 1
```

```
## Error in as.mcmc.list(x): object 'outlierTotal' not found
```

```r
gelman.diag(outlierTotal[3:4]) # the upper C.I is 1
```

```
## Error in as.mcmc.list(x): object 'outlierTotal' not found
```

## Manually manipulating MCMC chains

This section explains how to manually manipulate our MCMC chains.

### Discarding entire chains

Here we demonstrate how to discard chains that may not have converged.
Including chains that have not converged in downstream analysis is
likely to lead to nonsensical results. As an example, we demonstrate
how to remove the second chain from our domwnstream analysis.


```r
## Our chains look really good but let us discard chain 2, as an
## example in case we didn't believe it had converged.  It would be
## possible to remove more than one chain e.g. to remove 2 and 4 using
## c(2, 4).
removeChain <- 2
newTanMcmc <- tanTagm[seq_len(nChains)[-removeChain]]
```

```
## Error in tanTagm[seq_len(nChains)[-removeChain]]: object 'nChains' not found
```

```r
## Let check that it looks good
newTanMcmc
```

```
## Error in eval(expr, envir, enclos): object 'newTanMcmc' not found
```

```r
length(newTanMcmc) == (nChains - length(removeChain))
```

```
## Error in eval(expr, envir, enclos): object 'newTanMcmc' not found
```

```r
## Let have a look at our first chain
tanChain1 <- pRoloc:::chains(newTanMcmc)[[1]]
```

```
## Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]): there is no package called 'BiocInstaller'
```

```r
tanChain1@n ## Chain has 1500 iterations.
```

```
## Error in eval(expr, envir, enclos): object 'tanChain1' not found
```

We could repeat the convergence diagnostics of the previous section to
check for convergence again.

## Discarding iterations from each chain

We are happy that our chains have conveged. Let us recap where our
analysis up until this point. We have 3 chains each with 1500
iterations. We now demonstrate how to discard iterations from
individual chains, since they may have converged some number of
iteration into the analysis. Let us use only the iterations of last
half of the remaining three chains.  This also speeds up computations,
however too few iterations in the further analysis is likely to lead
to poor results. We find that using at least $1000$ iterations for
downstream analysis leads to stable results.



```r
## We need to clear this section up with a new function

n <- (tanChain1@n)/2 # Number of iterations to keep 750
```

```
## Error in eval(expr, envir, enclos): object 'tanChain1' not found
```

```r
K <- tanChain1@K # Number of components
```

```
## Error in eval(expr, envir, enclos): object 'tanChain1' not found
```

```r
N <- tanChain1@N # Number of Proteins
```

```
## Error in eval(expr, envir, enclos): object 'tanChain1' not found
```

```r
## Create storage for .MCMCChain
.MCMCChainlist <- vector("list", length = length(newTanMcmc))
```

```
## Error in vector("list", length = length(newTanMcmc)): object 'newTanMcmc' not found
```

```r
for(j in seq_len(length(newTanMcmc))) {

  tanChain <- pRoloc:::chains(newTanMcmc)[[j]]
  .ComponentParam <- tanChain@ComponentParam # This won't change

  ## Subset MCMC iterations
  retain <- seq.int(n + 1, tanChain@n) # retain 750 samples

  ## Check correct number of iterations
  stopifnot(ncol(tanChain@Component[, retain]) == n) # Second entry is 750

  ## Subset functions
  .Component <- tanChain@Component[, retain]
  .ComponentProb <- tanChain@ComponentProb[, retain, ]
  .Outlier <- tanChain@Outlier[, retain]
  .OutlierProb <- tanChain@OutlierProb[, retain, ]

  ## We can now create a new object of class MCMCChains
  ## make MCMCChain object
  .MCMCChainlist[[j]] <- pRoloc:::.MCMCChain(n = as.integer(n),
											 K = K,
											 N = N,
											 Component = .Component,
											 ComponentProb = .ComponentProb,
											 Outlier = .Outlier,
											 OutlierProb = .OutlierProb,
											 ComponentParam = .ComponentParam)

}
```

```
## Error in eval(expr, envir, enclos): object 'newTanMcmc' not found
```

```r
## Construct class MCMCChains
.ans <- pRoloc:::.MCMCChains(chains = .MCMCChainlist)
```

```
## Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]): there is no package called 'BiocInstaller'
```

```r
tanTagmparams <- pRoloc:::.MCMCParams(method = "TAGM.MCMC",
									  chains = .ans,
									  priors = tanTagm@priors,
									  summary = pRoloc:::.MCMCSummary())
```

```
## Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]): there is no package called 'BiocInstaller'
```

tanTagmParams is now an object of class `MCMCParams` with 3 chains
each with 750 iterations.



```r
## Check tanTagmParams object
pRoloc:::chains(tanTagmparams)[[1]]
```

```
## Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]): there is no package called 'BiocInstaller'
```

```r
pRoloc:::chains(tanTagmparams)[[2]]
```

```
## Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]): there is no package called 'BiocInstaller'
```

```r
pRoloc:::chains(tanTagmparams)[[3]]
```

```
## Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]): there is no package called 'BiocInstaller'
```

## Procesing and summarising MCMC results

### Populating the summary slot

The summary slot of the `tanTagmparams` is currently empty, we can now
populate the summary slot of `tanTagmparams` using the
`tagmMcmcProcess` function.


```r
## This will automatically pool chains to produce summary (easy to
## create single summaries by subsetting)
tanTagmparams <- tagmMcmcProcess(tanTagmparams)
```

```
## Error in tagmMcmcProcess(tanTagmparams): could not find function "tagmMcmcProcess"
```

```r
## Let look at this object
summary(tanTagmparams@summary@posteriorEstimates)
```

```
## Error in summary(tanTagmparams@summary@posteriorEstimates): object 'tanTagmparams' not found
```

For a sanity check, let us re-check the diagnostics.  This is
re-computed when we excute the `tagmMcmcProcess` function and located
in a `diagnostics` slot.



```r
## Recomputed diagnostics
tanTagmparams@summary@diagnostics
```

```
## Error in eval(expr, envir, enclos): object 'tanTagmparams' not found
```

Let us look at a summary of the analysis.



```r
summary(tanTagmparams@summary@tagm.joint)
```

```
## Error in summary(tanTagmparams@summary@tagm.joint): object 'tanTagmparams' not found
```

### Appending results to MSnSet

The `pRoloc` function `tagmPredict` can be used to append protein MCMC
results to the feature data of our object of class `MSnSet`. This
creates new columns in the feature data of the `MSnSet`, which can be
used for final analysis of our data.


```r
## We can now use tagmPredict
tan2009r1 <- tagmPredict(tan2009r1, params = tanTagmparams)
```

```
## Error in tagmPredict(tan2009r1, params = tanTagmparams): could not find function "tagmPredict"
```

## Visualising MCMC results

### Visualising prediction results

Now that we have processed our chains, checked convergence and
summarised the results into our `MSnset`; we can interrogate our data
for novel biological results. We use the `plot2D` function to view the
probabilitic allocations of proteins to sub-cellular niches.


```r
## Create prediction point size
ptsze <- exp(fData(tan2009r1)$tagm.mcmc.probability) - 1
```

```
## Error in fData(tan2009r1): object 'tan2009r1' not found
```

```r
## Create plot2D with pointer scaled with probability
plot2D(tan2009r1, fcol = "tagm.mcmc.allocation", cex = ptsze,
	   main = "protein pointer scaled with posterior localisation probability")
```

```
## Error in plot2D(tan2009r1, fcol = "tagm.mcmc.allocation", cex = ptsze, : could not find function "plot2D"
```

```r
addLegend(object = tan2009r1, where = "topleft", cex = 0.5)
```

```
## Error in addLegend(object = tan2009r1, where = "topleft", cex = 0.5): could not find function "addLegend"
```

### Visualising allocation uncertainty

By using the Shannon entropy we can globally visualise
uncertainty. Proteins can be scaled with their Shannon entropy and we
note that proteins with high Shannon entropy have high uncertainty.


```r
## Visualise shannon entropy
## Create prediction point size
ptsze2 <- 3 * fData(tan2009r1)$tagm.mcmc.mean.shannon
```

```
## Error in fData(tan2009r1): object 'tan2009r1' not found
```

```r
plot2D(tan2009r1, fcol = "tagm.mcmc.allocation", cex = ptsze2,
	   main = "protein pointer scaled with Shannon entropy")
```

```
## Error in plot2D(tan2009r1, fcol = "tagm.mcmc.allocation", cex = ptsze2, : could not find function "plot2D"
```

```r
addLegend(object = tan2009r1, where = "topleft", cex = 0.5)
```

```
## Error in addLegend(object = tan2009r1, where = "topleft", cex = 0.5): could not find function "addLegend"
```

### Extracting proteins of interest

Our data can be interrogated in other ways. We might be interested in
all the proteins that were confidently assigned as outliers. A GO
enrichment analysis could be performed on these proteins, since they
may reveal biologically interesting results.


```r
## Get outlier lists proteins with probability greater than 0.95 of being outlier
outliers <- rownames(tan2009r1)[fData(tan2009r1)$tagm.mcmc.outlier > 0.95]
```

```
## Error in rownames(tan2009r1): object 'tan2009r1' not found
```

```r
outliers
```

```
## Error in eval(expr, envir, enclos): object 'outliers' not found
```

### Extracting information for individual proteins

We might be interested in analysing the localisation of individual
proteins and interrogating which other potential localisations these
proteins might have.  A violin plot of the probabilistic allocation of
the following protein Q9VCK0 is visualised. This demonstrates and
quantifies the uncertainty in the allocation of this protein.


```r
## We can make this into a function
Q9VCK0 <- as.data.frame(tanChain@ComponentProb["Q9VCK0",,])
```

```
## Error in as.data.frame(tanChain@ComponentProb["Q9VCK0", , ]): object 'tanChain' not found
```

```r
colnames(Q9VCK0) <- getMarkerClasses(tan2009r1)
```

```
## Error in getMarkerClasses(tan2009r1): could not find function "getMarkerClasses"
```

```r
Q9VCK0melt <- melt(Q9VCK0)
```

```
## Error in melt(Q9VCK0): object 'Q9VCK0' not found
```

```r
colnames(Q9VCK0melt) <- c("Organelle","Probability")
```

```
## Error in colnames(Q9VCK0melt) <- c("Organelle", "Probability"): object 'Q9VCK0melt' not found
```

```r
gg2 <- ggplot(Q9VCK0melt,
			  aes(Organelle, Probability, width = (Probability))) +
	geom_violin(aes(fill = Organelle), scale = "width")
```

```
## Error in ggplot(Q9VCK0melt, aes(Organelle, Probability, width = (Probability))): object 'Q9VCK0melt' not found
```

```r
gg2 <- gg2 + theme_bw() +
	scale_fill_manual(values = getStockcol()[1:14]) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1),
		  axis.title.x = element_blank())
```

```
## Error in eval(expr, envir, enclos): object 'gg2' not found
```

```r
gg2 <- gg2 +
	ylab("Membership Probability") +
	ggtitle(paste0("Distribution of Subcellular Membership for Protein Q9VCK0" ))
```

```
## Error in eval(expr, envir, enclos): object 'gg2' not found
```

```r
gg2 <- gg2 +
	theme(legend.position="none")
```

```
## Error in eval(expr, envir, enclos): object 'gg2' not found
```

```r
print(gg2)
```

```
## Error in print(gg2): object 'gg2' not found
```
