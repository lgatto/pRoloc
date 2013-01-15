`pRoloc` is a Bioconductor package for the analysis of organelle proteomics data analysis.
It is available from Bioconductor `>= 2.12`.

The preferred installation procedure uses the Bioconductor infrastructure:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("pRoloc")
```

If you do not have a recent R version, you can install 
the `pRoloc` dependencies as follow and install `pRoloc` manually.

```r
deps <- c("MSnbase", "MLInterfaces", "mclust", "MSBVAR",
	      "caret", "e1071", "sampling", "class", "kernlab",
		  "nnet", "randomForest", "proxy", "BiocGenerics",
		  "RColorBrewer", "scales", "pRolocdata")
source("http://proteome.sysbiol.cam.ac.uk/lgatto/src/installPackages.R")
installPackages("pRoloc")
```

Download the appropriate package from the [Bioconductor landing page](http://www.bioconductor.org/packages/devel/bioc/html/pRoloc.html)
and install manually using `install.packages(..., repos = "NULL")` or the GUI front-end of your favourite R interface.
