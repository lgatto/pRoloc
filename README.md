[![Build Status](https://travis-ci.org/lgatto/pRoloc.svg?branch=master)](https://travis-ci.org/lgatto/pRoloc) [![codecov.io](https://codecov.io/github/lgatto/pRoloc/coverage.svg?branch=master)](https://codecov.io/github/lgatto/pRoloc?branch=master)

# A unifying bioinformatics framework for spatial proteomics

<img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/pRoloc/pRoloc.png" height="100"><img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/pRoloc/pRolocdata.png" height="100"><img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/pRoloc/pRolocGUI.png" height="100">

The `pRoloc` suite set of software are distributed as part of the
R/[Bioconductor](http://bioconductor.org/) project and are developed
at the [Computational Proteomics Unit](http://cpu.sysbiol.cam.ac.uk/)
and
[Cambridge Centre for Proteomics](http://proteomics.bio.cam.ac.uk/)
labs, at the University of Cambridge. See the Bioconductor landing
pages
([release](http://www.bioconductor.org/packages/release/bioc/html/pRoloc.html) and
[devel](http://www.bioconductor.org/packages/devel/bioc/html/pRoloc.html))
for the Bioc build status.



## Getting started

The `pRoloc` software comes with ample
documentation. The
[main tutorial](https://lgatto.github.io/pRoloc/articles/pRoloc-tutorial.html) provides
a broad overview of the package and its functionality.  See the
*Articles* tab for additional manuals.

Two associated
packages,
[`pRolocdata`](http://www.bioconductor.org/packages/devel/data/experiment/html/pRolocdata.html) and
[`pRolocGUI`](http://www.bioconductor.org/packages/release/bioc/html/pRolocGUI.html) offer
spatial proteomics data and a graphical user interface to
interactively explore the data.

Here are a set of
[video tutorial](https://www.youtube.com/playlist?list=PLvIXxpatSLA2loV5Srs2VBpJIYUlVJ4ow)
that illustrate the `pRoloc` framework.

## Help

Post your questions on the
[Bioconductor support site](https://support.bioconductor.org/),
tagging it with the package name `pRoloc` (the maintainer will
automatically be notified by email). If you identify a bug or have a
feature request, please open an
[issue](https://github.com/lgatto/pRoloc/issues) on the github
development page.

## Installation

The preferred installation procedure uses the Bioconductor
infrastructure:

```c
source("http://bioconductor.org/biocLite.R")
biocLite("pRoloc")
biocLite("pRolocdata")
biocLite("pRolocGUI")
```  

### Pre-release/development version

The pre-release/development code on github can be installed using
`biocLite`. Note that this requires a working R build environment (i.e
`Rtools` on Windows - see
[here](https://github.com/lgatto/teachingmaterial/wiki/R-package)). New
pre-release features might not be documented not thoroughly tested and
could substantially change prior to release. Use at your own risks.


```c
## install from github
biocLite("lgatto/pRoloc")
biocLite("lgatto/pRolocdata")
biocLite("ComputationalProteomicsUnit/pRolocGUI")
```

## References:

For refences about the software, how to use it and spatial proteomics
data analysis:

* Breckels LM, Mulvey CM, Lilley KS and Gatto L. A Bioconductor
  workflow for processing and analysing spatial proteomics
  data. F1000Research 2016,
  5:2926
  [doi:10.12688/f1000research.10411.1](https://f1000research.com/articles/5-2926/).

* Gatto L, Breckels LM, Burger T, Nightingale DJ, Groen AJ, Campbell
  C, Nikolovski N, Mulvey CM, Christoforou A, Ferro M, Lilley KS. A
  foundation for reliable spatial proteomics data analysis. Mol Cell
  Proteomics. 2014 Aug;13(8):1937-52. doi:
  10.1074/mcp.M113.036350. Epub 2014
  May 20. [PubMed PMID: 24846987](http://www.ncbi.nlm.nih.gov/pubmed/24846987)

* Gatto L, Breckels LM, Wieczorek S, Burger T, Lilley
  KS. Mass-spectrometry-based spatial proteomics data analysis using
  pRoloc and pRolocdata. Bioinformatics. 2014 May 1;30(9):1322-4. doi:
  10.1093/bioinformatics/btu013. Epub 2014
  Jan 11. [PubMed PMID: 24413670](http://www.ncbi.nlm.nih.gov/pubmed/24413670).

Specific algorithms available in the software:

* Breckels LM, Gatto L, Christoforou A, Groen AJ, Lilley KS, Trotter
  MW. The effect of organelle discovery upon sub-cellular protein
  localisation. J Proteomics. 2013 Aug 2;88:129-40. doi:
  10.1016/j.jprot.2013.02.019. Epub 2013
  Mar 21. [PubMed PMID: 23523639](http://www.ncbi.nlm.nih.gov/pubmed/23523639).

* Breckels LM, Holden S, Wojnar D, Mulvey CMM, Christoforou A, Groen
  AJ, Kohlbacher O, Lilley KS and Gatto L. Learning from heterogeneous
  data sources: an application in spatial proteomics. 2015 biorXiv,
  doi: http://dx.doi.org/10.1101/022152



#### More resource

* R and Bioconductor for proteomics
  [web page](http://lgatto.github.io/RforProteomics/) and
  [package](http://www.bioconductor.org/packages/release/data/experiment/html/RforProteomics.html)
  
* Bioconductor proteomics [workflow](http://bioconductor.org/help/workflows/proteomics/)

## Contributing

Contributions to the package are more than welcome. If you want to
contribute to this package, you should follow the same conventions as
the rest of the functions whenever it makes sense to do so. Feel free
to get in touch (preferable opening a
[github issue](https://github.com/lgatto/pRoloc/issues/)) to discuss
any suggestions. 

Please note that this project is released with a
[Contributor Code of Conduct](https://github.com/lgatto/pRoloc/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.
