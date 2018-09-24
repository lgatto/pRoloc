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
## Object of class "MCMCParams"
## Method: TAGM.MCMC
## Number of chains: 4
```

We now load the example data for which we performed this Bayesian
analysis. This Drosophila embryo spatial proteomics experiment is
from [@Tan2009].


```r
data(tan2009r1) ## get data from pRolocdata
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
nChains <- length(tanTagm)
nChains
```

```
## [1] 4
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
## iterations or entire chains: compute the number of outliers for
## each iteration for each chain
out <- mcmc_get_outliers(p_unst)
```

```
## Error in inherits(x, "MCMCParams"): object 'p_unst' not found
```

```r
## Using coda S3 objects to produce trace plots and histograms
plot(out[[1]], col = "blue", main = "Chain 1")
```

```
## Error in plot(out[[1]], col = "blue", main = "Chain 1"): object 'out' not found
```

```r
plot(out[[2]], col = "red", main = "Chain 2")
```

```
## Error in plot(out[[2]], col = "red", main = "Chain 2"): object 'out' not found
```

```r
plot(out[[3]], col = "green", main = "Chain 3")
```

```
## Error in plot(out[[3]], col = "green", main = "Chain 3"): object 'out' not found
```

```r
plot(out[[4]], col = "orange", main = "Chain 4")
```

```
## Error in plot(out[[4]], col = "orange", main = "Chain 4"): object 'out' not found
```

We can use the *[coda](https://CRAN.R-project.org/package=coda)* package to produce
summaries of our chains. Here is the `coda` summary for the first
chain.


```r
## all chains average around 360 outliers
summary(out[[1]])
```

```
## Error in summary(out[[1]]): object 'out' not found
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
gelman.diag(out) ## the Upper C.I. is 1 so mcmc has clearly converged
```

```
## Error in as.mcmc.list(x): object 'out' not found
```

We can also look at the Gelman diagnostics statistics for pairs of chains.


```r
## We can also check individual pairs of chains for convergence
gelman.diag(out[1:2]) # the upper C.I is 1.01
```

```
## Error in as.mcmc.list(x): object 'out' not found
```

```r
gelman.diag(out[c(1,3)]) # the upper C.I is 1
```

```
## Error in as.mcmc.list(x): object 'out' not found
```

```r
gelman.diag(out[3:4]) # the upper C.I is 1
```

```
## Error in as.mcmc.list(x): object 'out' not found
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
newTanMcmc <- tanTagm[-2]

## Let check that it looks good
newTanMcmc
```

```
## Object of class "MCMCParams"
## Method: TAGM.MCMC
## Number of chains: 3
```

```r
length(newTanMcmc) == (nChains - 1)
```

```
## [1] TRUE
```

```r
## Let have a look at our first chain
tanChain1 <- pRoloc:::chains(newTanMcmc)[[1]]
tanChain1@n ## Chain has 1500 iterations.
```

```
## [1] 1500
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
n <- (tanChain1@n)/2 # Number of iterations to keep 750
tanTagm2 <- mcmc_burn_chains(newTanMcmc, n)
```

`tanTagm2` is now an object of class `MCMCParams` with 3 chains each
with each 750 iterations.



```r
## Check tanTagmParams object
pRoloc:::chains(tanTagm2)[[1]]
```

```
## Object of class "MCMCChain"
##  Number of components: 11
##  Number of proteins: 677
##  Number of iterations: 750
```

```r
pRoloc:::chains(tanTagm2)[[2]]
```

```
## Object of class "MCMCChain"
##  Number of components: 11
##  Number of proteins: 677
##  Number of iterations: 750
```

```r
pRoloc:::chains(tanTagm2)[[3]]
```

```
## Object of class "MCMCChain"
##  Number of components: 11
##  Number of proteins: 677
##  Number of iterations: 750
```

## Processing and summarising MCMC results

### Populating the summary slot

The summary slot of the `tanTagm2` is currently empty, we can now
populate the summary slot of `tanTagm2` using the `tagmMcmcProcess`
function.


```r
## This will automatically pool chains to produce summary (easy to
## create single summaries by subsetting)
tanTagm2 <- tagmMcmcProcess(tanTagm2)

## Let look at this object
summary(tanTagm2@summary@posteriorEstimates)
```

```
##       tagm.allocation tagm.probability tagm.probability.notOutlier
##  ER           :196    Min.   :0.3055   Min.   :0.0002031
##  PM           :192    1st Qu.:0.8439   1st Qu.:0.1712611
##  Ribosome 40S : 75    Median :0.9810   Median :0.6001312
##  mitochondrion: 61    Mean   :0.8927   Mean   :0.5278532
##  Proteasome   : 49    3rd Qu.:0.9980   3rd Qu.:0.8603698
##  Ribosome 60S : 32    Max.   :1.0000   Max.   :0.9982244
##  (Other)      : 72
##  tagm.probability.Outlier tagm.probability.lowerquantile
##  Min.   :0.001776         Min.   :0.001053
##  1st Qu.:0.139630         1st Qu.:0.591425
##  Median :0.399869         Median :0.942714
##  Mean   :0.472147         Mean   :0.769211
##  3rd Qu.:0.828739         3rd Qu.:0.995483
##  Max.   :0.999797         Max.   :0.999988
##
##  tagm.probability.upperquantile tagm.mean.shannon
##  Min.   :0.4961                 Min.   :0.0000254
##  1st Qu.:0.9698                 1st Qu.:0.0141329
##  Median :0.9959                 Median :0.0911157
##  Mean   :0.9655                 Mean   :0.2386502
##  3rd Qu.:0.9994                 3rd Qu.:0.4050057
##  Max.   :1.0000                 Max.   :1.3051988
##
```

For a sanity check, let us re-check the diagnostics.  This is
re-computed when we excute the `tagmMcmcProcess` function and located
in a `diagnostics` slot.



```r
## Recomputed diagnostics
tanTagm2@summary@diagnostics
```

```
##              Point est. Upper C.I.
## outlierTotal    1.00021   1.000274
```

Let us look at a summary of the analysis.



```r
summary(tanTagm2@summary@tagm.joint)
```

```
##   Cytoskeleton             ER             Golgi
##  Min.   :0.0000000   Min.   :0.0000   Min.   :0.0000000
##  1st Qu.:0.0000000   1st Qu.:0.0000   1st Qu.:0.0000000
##  Median :0.0000040   Median :0.0000   Median :0.0000002
##  Mean   :0.0227762   Mean   :0.2826   Mean   :0.0503043
##  3rd Qu.:0.0004646   3rd Qu.:0.9232   3rd Qu.:0.0046922
##  Max.   :0.7798512   Max.   :0.9993   Max.   :0.9993792
##     Lysosome         mitochondrion          Nucleus
##  Min.   :0.0000000   Min.   :0.0000000   Min.   :0.0000000
##  1st Qu.:0.0000001   1st Qu.:0.0000000   1st Qu.:0.0000000
##  Median :0.0000016   Median :0.0000000   Median :0.0000000
##  Mean   :0.0167424   Mean   :0.0916050   Mean   :0.0365126
##  3rd Qu.:0.0000793   3rd Qu.:0.0000158   3rd Qu.:0.0000059
##  Max.   :0.9795195   Max.   :0.9998841   Max.   :0.9999886
##    Peroxisome              PM             Proteasome
##  Min.   :0.0000000   Min.   :0.000000   Min.   :0.0000000
##  1st Qu.:0.0000000   1st Qu.:0.000000   1st Qu.:0.0000000
##  Median :0.0000002   Median :0.000003   Median :0.0000002
##  Mean   :0.0055262   Mean   :0.278256   Mean   :0.0726235
##  3rd Qu.:0.0000232   3rd Qu.:0.764944   3rd Qu.:0.0005521
##  Max.   :0.8831085   Max.   :0.999998   Max.   :0.9979738
##   Ribosome 40S       Ribosome 60S
##  Min.   :0.000000   Min.   :0.0000000
##  1st Qu.:0.000000   1st Qu.:0.0000000
##  Median :0.000000   Median :0.0000001
##  Mean   :0.096071   Mean   :0.0470258
##  3rd Qu.:0.001051   3rd Qu.:0.0001129
##  Max.   :0.996639   Max.   :0.9885391
```

### Appending results to an MSnSet

The `pRoloc` function `tagmPredict` can be used to append protein MCMC
results to the feature data of our object of class `MSnSet`. This
creates new columns in the feature data of the `MSnSet`, which can be
used for final analysis of our data.


```r
## We can now use tagmPredict
tan2009r1 <- tagmPredict(tan2009r1, params = tanTagm2)
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

## Create plot2D with pointer scaled with probability
plot2D(tan2009r1,
       fcol = "tagm.mcmc.allocation",
       cex = ptsze,
       main = "protein pointer scaled with posterior localisation probability")

addLegend(object = tan2009r1, where = "topleft", cex = 0.5)
```

![plot of chunk mcmc-vis](figure/mcmc-vis-1.png)

### Visualising allocation uncertainty

By using the Shannon entropy we can globally visualise
uncertainty. Proteins can be scaled with their Shannon entropy and we
note that proteins with high Shannon entropy have high uncertainty.


```r
## Visualise shannon entropy
## Create prediction point size
ptsze2 <- 3 * fData(tan2009r1)$tagm.mcmc.mean.shannon
plot2D(tan2009r1, fcol = "tagm.mcmc.allocation", cex = ptsze2,
       main = "protein pointer scaled with Shannon entropy")
addLegend(object = tan2009r1, where = "topleft", cex = 0.5)
```

![plot of chunk mcmc-vis2](figure/mcmc-vis2-1.png)

### Extracting proteins of interest

Our data can be interrogated in other ways. We might be interested in
all the proteins that were confidently assigned as outliers. A GO
enrichment analysis could be performed on these proteins, since they
may reveal biologically interesting results.


```r
## Get outlier lists proteins with probability greater than 0.95 of being outlier
outliers <- rownames(tan2009r1)[fData(tan2009r1)$tagm.mcmc.outlier > 0.95]
outliers
```

```
##  [1] "P04412"   "Q7KJ73"   "Q00174"   "Q9Y105"   "Q9VTZ5"   "B7Z0C1"
##  [7] "P46150"   "Q9VN14"   "Q9VU58"   "Q9V498"   "Q9V496"   "P06607"
## [13] "Q7K0F7"   "Q9V8R9"   "Q9VT75"   "A8JNJ6"   "Q8SZ38"   "Q960V7"
## [19] "Q7KA66"   "P98163"   "Q9VE24"   "Q9VE79"   "A8JV22"   "O96824"
## [25] "P23226"   "Q9VKF0"   "Q9VZL6"   "P11584"   "Q9VDK9"   "Q9V4A7"
## [31] "Q9W303"   "Q8MLV1"   "Q07152"   "Q9VEE9"   "Q9VXE5"   "P13395"
## [37] "Q7JYX2"   "P91938"   "Q9U969"   "Q7JZY1"   "O18388"   "P02844"
## [43] "NO_ID_7"  "P27716"   "Q9W478"   "Q8IN49"   "O18337"   "A1Z6L9"
## [49] "Q86P66"   "Q9VK04"   "Q9GU68"   "Q8INH5"   "O18333"   "P54351"
## [55] "P15215"   "P33450"   "Q24298"   "Q9W404"   "Q9VPQ2"   "Q9VF20"
## [61] "Q8IPU3"   "Q7JVX3"   "Q9VV60"   "Q9V427"   "Q7JQR3"   "P52295"
## [67] "A8JRB8"   "Q9VTU4"   "Q24007"   "B8A403"   "Q00963"   "Q94887"
## [73] "Q4AB31"   "Q95SY7"   "O15943"   "Q9VDV3"   "B7Z0D3"   "Q9VLT3"
## [79] "Q9VUQ7"   "Q9VAG2"   "A1Z6H7"   "Q27415"   "Q9VHL2"   "NO_ID_16"
## [85] "Q94518"   "Q7KN94"   "P29742"   "Q9VXB0"   "B7Z036"   "Q8IRH6"
## [91] "P17917"   "M9PDB2"   "Q9VVA0"   "Q9VHP0"
```

### Extracting information for individual proteins

We might be interested in analysing the localisation of individual
proteins and interrogating which other potential localisations these
proteins might have.  A violin plot of the probabilistic allocation of
the following protein Q9VCK0 is visualised. This demonstrates and
quantifies the uncertainty in the allocation of this protein.


```r
mcmc_plot_probs(tanTagm, "Q9VCK0", tan2009r1)
```

```
## Loading required package: reshape2
```

```
##
## Attaching package: 'reshape2'
```

```
## The following objects are masked from 'package:reshape':
##
##     colsplit, melt, recast
```

```
## No id variables; using all as measure variables
```

![plot of chunk mcmc-gg2](figure/mcmc-gg2-1.png)
