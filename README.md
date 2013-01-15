`pRoloc` is a Bioconductor package for the analysis of organelle proteomics data analysis.
It is available from Bioconductor >= 2.12.

The preferred installation procedure uses the Bioconductor infrastructure:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("pRoloc")
```

If you do not have a recent R version, you can install 
the `pRoloc` dependencies as follow and install `pRoloc` manually.

```r
source("http://proteome.sysbiol.cam.ac.uk/lgatto/src/getDependencies.R")
source("http://proteome.sysbiol.cam.ac.uk/lgatto/src/installDependencie.R")
installDependencies("pRoloc")
```

Download the appropriate package from the [Bioconductor landing page](http://www.bioconductor.org/packages/devel/bioc/html/pRoloc.html)
and install manually using `install.packages(..., repos = "NULL")` or other the GUI front-end of your favourite R interface.
