`pRoloc` is a Bioconductor package for the analysis of organelle proteomics data analysis.
It is available from Bioconductor `>= 2.12`.

The preferred installation procedure uses the Bioconductor infrastructure:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("pRoloc")
```  

If you do not have a recent R version (R-devel), you can install 
the `pRoloc` dependencies as follow and install `pRoloc` manually.

```r
deps <- c("MSnbase", "MLInterfaces", "mclust", "MSBVAR", 
          "caret", "e1071", "sampling", "class", "kernlab",
          "nnet", "randomForest", "proxy", "BiocGenerics",
          "RColorBrewer", "scales")
source("http://proteome.sysbiol.cam.ac.uk/lgatto/src/installPackages.R")
installPackages(deps)
```

Note that you will also need the devel version of `MSnbase`, available on its 
[Bioconductor landing page](http://www.bioconductor.org/packages/2.12/bioc/html/MSnbase.html).

Download the appropriate package from the [Bioconductor landing page](http://www.bioconductor.org/packages/devel/bioc/html/pRoloc.html)
and install manually using `install.packages(..., repos = "NULL")` or the GUI front-end of your favourite R interface.

Alternatively, to install the non-official git version, you can proceed as follows. This method requires the package building mechanism to be readily available on your system.

```r
library("devtools")
install_github("pRoloc", "lgatto")
```

You can also install the associated data package [`pRolocdata`](http://bioconductor.org/packages/devel/data/experiment/html/pRolocdata.html) 
with multiple test data sets.
