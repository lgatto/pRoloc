# pRoloc 1.47

## Changes in version 1.47.1

- Lisa Breckels is now maintainer.

## Changes in version 1.47.0

- New devel version

# pRoloc 1.45

## Changes in version 1.45.2

- `pRolocmarkers()` has a new `version` argument, to allow for new
  markers versions to be added.
- 14 new marker sets have been added to `pRolocmarkers()` under
  `version = 2` (new default).
- Documentation for `pRolocmarkers()` has been updated to include a
  description of `version = 2` markers.

## Changes in version 1.45.1

- Import 'mclust::estep*()'.

## Changes in version 1.45.0

- New devel version

# pRoloc 1.43

## Changes in version 1.43.2

- Fix/update dunkley2006param object.

## Changes in version 1.43.1

- Fix syntax in man pages.

## Changes in version 1.43.0

- New devel version

# pRoloc 1.39

## Changes in version 1.39.1
- Update transfer learning vignette to use hpar 1.41.

## Changes in version 1.39.0
- New devel version (Bioc 3.17)

# pRoloc 1.37

## Changes in version 1.37.1
- Fix bug in PerTubro classifiction function (see #146 and #147),
  contributed by mgerault.

# pRoloc 1.33

# pRoloc 1.31

## Changes in version 1.31.1
- Fix failing unit test, by setting RNGseed in `SerialParam()` (fix by
  ococrook, see #142).

## Changes in version 1.31.0
- New devel version (Bioc 3.14)

## Changes in version 1.31.3
- lopims() function moved to lgatto/lopims package on GitHub

## Changes in version 1.31.2
- Update dunkley2006params, andy2011params and MartInterfaces objects

## Changes in version 1.31.1
- Suggest magick (needed to build vignette)

## Changes in version 1.31.0
- New devel version (Bioc 3.13)

# pRoloc 1.29

## Changes in version 1.29.0
- New devel version (Bioc 3.12)

## Changes in version 1.29.1
- Removed duplicate colour in getStockCol and updated the colour
  palette <Tue 2020-09-08>

# pRoloc 1.28

## Changes in version 1.28.0
- New release version (Bioc 3.11)

# pRoloc 1.27

## Changes in version 1.27.6
- Depend on MLInterfaces 1.67.10

## Changes in version 1.27.5
- import stats4::plot <2020-03-26 Thu>

## Changes in version 1.27.4
- Import missing `mclust::mclustBIC`

## Changes in version 1.27.3
- Remove exprsToRatio man page (function in MSnbase and is deprecated)

## Changes in version 1.27.2
- Fix errors related to R-devel

## Changes in version 1.27.1
- Merged plotting payes PR

## Changes in version 1.27.0
- Version bump for Bioc 3.11 (devel)


# pRoloc 1.25

## Changes in version 1.25.3
- New spatial2D function <2019-09-24 Tue>

## Changes in version 1.25.2
- Fix new biomart attribute <2019-08-09 Fri>
- Bug fix: pass fcol to helper function (see
  https://support.bioconductor.org/p/123614/) <2019-08-09 Fri>

## Changes in version 1.25.1
- Always use `mvtnorm::dmvnorm` <2019-06-21 Fri>

## Changes in version 1.25.0
- Version bump for Bioc 3.10 (devel)

# pRoloc 1.24

## Changes in version 1.24.0
- Version bump for Bioc 3.9 (release)

# pRoloc 1.23

## Changes in version 1.23.4
- Remove plain NEWS now that R supports NEWS.md

## Changes in version 1.23.3
- Link to TAGM F1000 workflow <2019-04-11 Thu>
- Add ref to TAGM F1000 in CITATION and README <2019-04-11 Thu>
- Update biomart interface data and AnnotationParams-related code
  <2019-04-12 Fri>

## Changes in version 1.23.2
- Add mcmc helper functions.
- New logPosteriors accessor for MAPParams object.
- New plotConsProfiles function, contributed by Tom Smith (see PR
  #131).
- Update TAGM vignette, refering to workflow <2019-03-14 Thu>

## Changes in version 1.23.1
- Add note about missing data in QSep man <2018-10-31 Wed>
- Add mcmc-helper documentation <2018-11-04 Sun>
- Fixing numerical instability in Cholesky decomposition (see #124)
  <2018-11-30 Fri>
- Add option to display or not display grid in plot2D <2018-12-12 Wed>
- `mrkHClust` now uses `mrkConsProfiles` and returns the `hclust`
  object - see issue #130 for details and background <2018-12-19 Wed>

## Changes in version 1.23.0
- New version for Bioc 3.9 devel

# pRoloc 1.22

## Changes in version 1.22.0
- New version for Bioc 3.8 release

# pRoloc 1.21

## Changes in version 1.21.9
- Fix type in vignette header <2018-09-18 Tue>
- Fix bug in plot method for ThetaRegRes object <2018-09-24 Mon>

## Changes in version 1.21.8
- Add an `fcol` argument to `plotDist` to plot and colour all profiles
  <2018-08-09 Thu>

## Changes in version 1.21.7
- Use BiocManager in vignette
- Fix bug in plot2D: pass ... to hexbin <2018-08-02 Thu>

## Changes in version 1.21.6
- Use BiocManager in installation instructions

## Changes in version 1.21.5
- Added new section in Bayesian spatial proteomics vignette detailing
  mcmc output processing <2018-07-07 Sat>

## Changes in version 1.21.4
- Fix bugs in tagmMcmcPredict, where fcol was ignored <2018-06-05 Tue>
- Order vignettes by prefixing the files with numbers <2018-06-05 Tue>

## Changes in version 1.21.3
- New TAGM-MCMC generative model, contributed by Oliver Crook
  <2018-05-18 Fri>

## Changes in version 1.21.2
- Version bump for BiocStyle update: Vignette needed to be rebuilt to
  have bug fixed in BiocStyle footnote rendering.

## Changes in version 1.21.1
- Fix bug in higlightOnPlot with missing fcol (see #105)
  <2018-05-03 Thu>
- New TAGM-MAP generative model, contributed by Oliver Crook
  <2018-05-18 Fri>
- New `plotEllipse` function to visualise and assess TAGM models
  <2018-05-18 Fri>

# pRoloc 1.20

## Changes in version 1.20.0
- New Bioconductor release 3.7

# pRoloc 1.19

## Changes in version 1.19.4
- Fix regression bug in knntl function <2018-04-12 Thu>

## Changes in version 1.19.3
- Use `dplyr::left_join` without attaching `dplyr` to avoid collision
  between `Biobase::exprs` and `dplyr::exprs` <2018-04-04 Wed>.
- Typo in warning to install rgl <2018-03-27 Tue>

## Changes in version 1.19.2
- Fix bug in QSep that prevented to set non-default fcol
  <2018-01-29 Mon>

## Changes in version 1.19.1
- Fix bug in private dimred and set appropriate number of colnames
  <2017-11-07 Tue>
- New `nipals` method in dimensionality reduction for plot2D (closes
  issue #103) <2018-01-16 Tue>

## Changes in version 1.19.0
- Bioconductor devel 3.7

# pRoloc 1.18

## Changes in version 1.18.0
- Bioconductor release 3.6

# pRoloc 1.17

## Changes in version 1.17.5
- Filtering for unique features when running plot2D with t-SNE method
  <2017-10-15 Sun>

## Changes in version 1.17.4

- Added new (private) dimred function that computes dimensionality
  reduction <2017-06-05 Mon>
- Add F1000research workflow to citations <2017-06-22 Thu>
- Classification functions now return the classification score matrix
  for all classes as a single column in fData, rather that each class
  as its own fData column. <2017-09-01 Fri>

## Changes in version 1.17.3
- Convert vignettes to Rmd with html output <2017-05-25 Thu>
- Import, rather than suggest Rtsne <2017-05-25 Thu>

## Changes in version 1.17.2
- phenoDisco speed improvements and added support for t-SNE
  <2017-05-19 Fri>

## Changes in version 1.17.1
- Support Rtsne's new pca_center and pca_scale arguments
  <2017-05-02 Tue>

## Changes in version 1.17.0
- Version bump for Bioc devel 3.6

# MSnbase 1.16

## Changes in version 1.16.0
- Version bump for Bioc release 3.5

# MSnbase 1.15

## Changes in version 1.15.9
- update biomart attribute names to relect changes <2017-04-20 Thu>

## Changes in version 1.15.8
- New mrkConsProfiles function to calculate average/consensus marker
  profiles <2017-04-11 Tue>

## Changes in version 1.15.7

 - Fix warnings and notes <2017-02-25 Sat>
 - Add section about dimensionality methods reduction and t-SNE in the
   tutorial <2017-03-07 Tue>
 - Fix error due to new uniprot attribute names <2017-04-06 Thu>

## Changes in version 1.15.6

 - fix (unexported) remap function - see issue #92 <2016-12-15 Thu>
 - plot2Ds now only adds segments when the featureNames are identical
   <2017-01-11 Wed>
 - Import (rather than Suggest) hexbin <2017-01-18 Wed>
 - Increase margin in QSep plotting (contributed by S. Gibb)
   <2017-02-07 Tue>

## Changes in version 1.15.5

 - Update plotDist to use sampleNames() to label x axis ticks and
   getUnknowncol() as default pcol (see issue #91) <2016-11-08 Tue>
 - xlab and ylab are args in plotDist <2016-12-04 Sun>

## Changes in version 1.15.4

 - Update human markers - see pRolocdata's issue 21
   https://github.com/lgatto/pRolocdata/issues/21 <2016-11-08 Tue>

## Changes in version 1.15.3

 - Fix bug in plot2D to ignore fcol when using hexbin method
        <2016-11-04 Fri>

## Changes in version 1.15.2

 - Fix Arabidopsis parameters for biomaRt questions: 'TAIR locus ID'
   is now 'Stable gene ID'. Changed in several manual files and
   dunkley2006params. <2016-11-02 Wed>

## Changes in version 1.15.1

 - new plot3D function <2016-10-27 Thu>
 - update CITATION <2016-10-27 Thu>
 - use predict:::predict.plsa, as it is not exported anymore
   <2016-11-01 Tue>

# MSnbase 1.13

## Changes in version 1.13.16

 - removing visualTest from suggest, as plotting tests currently
   disabled (see https://github.com/MangoTheCat/visualTest/issues/19 for
   details) <2016-10-07 Fri>

## Changes in version 1.13.15

 - Fix bug in plot,QSep with NAs <2016-09-09 Fri>
 - Adapt tl vignette to new hpar <2016-09-15 Thu>

## Changes in version 1.13.14

 - Pass ... in levelPlot,QSep <2016-09-06 Tue>
 - Various improvements to the tutorial vignette <2016-09-07 Wed>
 - Fix bug in QSep to honour fcol <2016-09-09 Fri>

## Changes in version 1.13.13

 - Added plot3D, equivalent to plot2D, but using rgl in 3 dimensions
   <2016-09-03 Sat>
 - Parametrise logging in knntlOptimisation to reduce object size
   <2016-09-06 Tue>
 - New mrkHClust function to plot a dendrogram of subcellular markers
   <2016-09-06 Tue>

## Changes in version 1.13.12

 - Remove verbose argument from getMarkerClasses function
   <2016-08-23 Tue>
 - plot2D has new addLegend argument <2016-08-23 Tue>
 - addLegend default position is now bottomleft (requested by Lisa)
   <2016-08-23 Tue>
 - new hexbin plot2D method <2016-08-28 Sun>

## Changes in version 1.13.11

 - Add a norm agument to QSep's plotting functions to visualise
   normalised and raw distances <2016-08-22 Mon>
 - New zeroInBinMSnSet function to visualise rowSums <2016-08-23 Tue>

## Changes in version 1.13.10

 - Send biomart queries in getGOFromFeatues in chunks - see issue #85
   for details <2016-08-20 Sat>

## Changes in version 1.13.9

 - new QSep infrastructure to assess spatial resolution
   <2016-08-17 Wed>
 - Use %age of total variance in plot2D's scree plot <2016-08-17 Wed>
 - fixed bug in combineThetaRegRes <2016-08-19 Fri>

## Changes in version 1.13.8

 - fDataToUnknown accepts from = NA <2016-08-10 Wed>

## Changes in version 1.13.7

 - Added an lda method to plot2D <2016-08-08 Mon>

## Changes in version 1.13.6

 - Fixed bug in addGoAnnotations <2016-07-29 Fri>

## Changes in version 1.13.5

 - Fixed bug in show method for ThetaRegRes <2016-06-08 Wed>
 - Profiled knntl code so it's now much faster. Added more unit tests
   for knntl <2016-06-10 Fri>
 - Fixed bug in plotDist x-axis labelling <2016-06-14 Tue>

## Changes in version 1.13.4

 - Updated mouse pRolocmarkers <2016-06-02 Thu>

## Changes in version 1.13.3

 - Fix bug in mirrorX/Y highlightOnPlot - see problem 1 in issue #79 <2016-05-31 Tue>

## Changes in version 1.13.2

 - Version bump to trigger package rebuilding now that purl()'ing
   issue has been correctly identified. knitr does not create
   purl()'ed (Stangle equivalent) .R files if _R_CHECK_TIMINGS_ is
   set, which the build system was setting. Now it's not set, so these
   .R files are now created. See
   https://github.com/yihui/knitr/issues/1212 for more. d.tenenbaum

## Changes in version 1.13.1

 - Added keepNA argument to goTermToId so that if a GO term becomes
   obsolete and you cannot replace it with the ID name, you have
   the option to either replace it with a NA (previous and current
   default option) or now with keepNA = FALSE the term name will be
   replaced with the id name <2016-05-25 Wed>

 - Bump version of all packages that use knitr for vignettes. This is
   because of an issue (now fixed) in knitr which failed to create
   purl()'ed R files from vignette sources and include them in the
   package. This version bump will cause these packages to propagate
   with those R files included. d.tenenbaum


## Changes in version 1.13.0

 - Version bump for new Bioc devel

# MSnbase 1.12

## Changes in version 1.12.0

 - Version bump for Bioc release version 3.3

# MSnbase 1.11

## Changes in version 1.11.23

 - New gomarkers functionality for adding annotation information
   to spatial proteomics data and accompanying new vignette
   <2016-04-18 Mon>
 - Added unit tests <2016-04-20 Wed> <2016-04-21 Thu>
 - Moved makeNaData[2] and whichNA to MSnbase <2016-04-21 Thu>
 - Renamed addGoMarkers and orderGoMarkers to addGoAnnotations
   and orderGoAnnotations and all associated documentation.
   Vignette also renamed to pRoloc-goannotations
   <2016-04-21 Thu>

## Changes in version 1.11.22

 - fix bug in knntlOptimisation to allow passing of th matrix
   with 1 column <2016-04-08 Fri>

## Changes in version 1.11.21

 - fix bug in getParams method = 'max' <2016-04-07 Thu>

## Changes in version 1.11.20

 - Update plotDist signature to support different types and pch
   <2016-04-05 Tue>

## Changes in version 1.11.19

  - Update dunkley2006params <2016-04-01 Fri>

## Changes in version 1.11.18

 - Update dunkley2006params <2016-03-30 Wed>

## Changes in version 1.11.17

 - Selective imports <2016-03-20 Sun>

## Changes in version 1.11.19

  - Update dunkley2006params <2016-04-01 Fri>

## Changes in version 1.11.18

 - Update dunkley2006params <2016-03-30 Wed>

## Changes in version 1.11.17

 - Selective imports <2016-03-20 Sun>

## Changes in version 1.11.16

 - Update package startup msg <2016-03-11 Fri>

## Changes in version 1.11.15

 - Clarify that score during optim are not a reflection of final
   assignment accuracy <2016-03-01 Tue>

## Changes in version 1.11.14

 - fix build error due to doi/url confusion <2016-03-01 Tue>

## Changes in version 1.11.13

 - new method argument added to knntlOptimisation that allows
   optimisation of class weights as per Wu and Dietterich's
   original k-NN TL method <2016-02-19 Fri>
 - seed argument added to knntlOptimisation for reproducibility
   <2016-02-22 Mon>
 - New section in tl vignette describing preparation of auxiliary PPI
   data <2016-02-29 Mon>

## Changes in version 1.11.12

 - the colour and pch setters now invisibly return the old values
   <2016-02-17 Wed>

## Changes in version 1.11.11

 - added error message when sampleNames differ between datasets
   in an MSnSetList when using remap function <2016-02-11 Thu>

## Changes in version 1.11.10

 - Update colours man page to document change in default colour
   palette <2016-02-08 Mon>

## Changes in version 1.11.9

 - Lisa's colour palette is now default. Old colours can be accessed
   and set with get/setOldcol <2016-02-04 Thu>

## Changes in version 1.11.8

 - New Lisa cols and changed default unknown col <2016-02-03 Wed>
 - mrkVecToMat has been updated so that the column order reflects
   the factor levels of fcol, rather than calling unique on fcol.
   This change means that the order of the classes in fcol are
   now consistent between plot2D and new visualisation apps that
   rely on mrkVecToMat. <2016-02-03 Wed>

## Changes in version 1.11.7

 - Various non-visible changes. <2016-01-19 Tue>

## Changes in version 1.11.6

 - new classWeights function <2015-12-26 Sat>

## Changes in version 1.11.5

 - highlightOnPlot support labels = TRUE to use featureNames as labels
   <2015-12-21 Mon>
 - selective ggplot2 import <2015-12-21 Mon>
 - highlightOnPlot also support a vector of feature names in addition
   to an instance of class FeaturesOfInterest <2015-12-21 Mon>

## Changes in version 1.11.4

 - plot2D: Mirror PCs even when not plotting, addressing issue #68
   <2015-12-18 Fri>

## Changes in version 1.11.3

 - Update dunkley2006params to use plant_mart_30 <2015-12-16 Wed>
 - API change in plot2D: to plot data as is, i.e. without any
   transformation, method can be set to "none" (as opposed to passing
   pre-computed values to method as in previous versions). If object
   is an MSnSet, the untransformed values in the assay data will be
   plotted. If object is a matrix with coordinates, then a matching
   MSnSet must be passed to methargs. <2015-12-16 Wed>

## Changes in version 1.11.2

 - Internally using MartInterface to query individual mart servers to
   bypass the biomart.org downtime <2015-12-09 Wed>

## Changes in version 1.11.1

 - New orgQuants function and update to getPredictions
   <2015-10-13 Tue>
 - Deprecate minClassScore replaced by getPredictions
   <2015-10-19 Mon>
 - Add pRolocVisMethods and check for new apps in
   pRolocGUI <2015-10-19 Mon>
 - new fDataToUnknown function <2015-10-23 Fri>
 - New section in vignette describing readMSnSet2 <2015-11-30 Mon>

## Changes in version 1.11.0

 - Bioc devel version 3.3

# MSnbase 1.10

## Changes in version 1.10.0

 - Bioc release version 3.2

# MSnbase 1.9

## Changes in version 1.9.7

 - New SpatProtVis visualisation class <2015-08-13 Thu>
 - add link to explanation of supportive/uncertain reliability scores
   in tl vignette <2015-09-02 Wed>

## Changes in version 1.9.6

 - Update REAMDE with TL ref
 - Update refs in lopims documentation <2015-07-30 Thu>

## Changes in version 1.9.5

 - update doc <2015-07-15 Wed>

## Changes in version 1.9.4

 - Add reference to TL paper and link to lpSVM code <2015-07-06 Mon>
 - highlightOnPlot throws a warning and invisibly returns NULL instead
   of an error when no features are in the object <2015-07-08 Wed>
 - highlightOnPlot has a new labels argument <2015-07-10 Fri>

## Changes in version 1.9.3

 - Clarify error when no annotation params are provided
   <2015-05-11 Mon>
 - support for matrix-encoded markers <2015-05-19 Tue>
 - New default in addLegend: bty = "n" <2015-05-20 Wed>
 - getMarkers now supports matrix markers <2015-05-20 Wed>
 - getMarkerClasses now supports matrix markers <2015-05-20 Wed>
 - markerMSnSet and unknownMSnSet now support matrix markers
   <2015-05-20 Wed>
 - sampleMSnSet now supports matrix markers <2015-05-23 Sat>
 - updated yeast markers and added uniprot ids <2015-05-27 Wed>
 - plot2D support a pre-calculated dim-reduced data matrix as method
   parameter to avoid recalculation <2015-05-27 Wed>

## Changes in version 1.9.2

 - Lisa's colour palette <2015-05-08 Fri>

## Changes in version 1.9.1

 - new plot2Ds function to overlay two data sets on the same PCA plot
   [2015-04-17 Fri]
 - regenerate biomart data used by setAnnotationParams
   [2015-04-24 Fri]
 - new setStockcolGui function to set the default colours manually via
   a simple interface [2015-04-29 Wed]
 - new move2Ds function to produce an transition movie between two
   MSnSets [2015-04-29 Wed]
 - functions to convert GO ids to/from terms. See ?goTermToId for
   details <2015-05-08 Fri>

## Changes in version 1.9.0

 - new devel version for Bioc 3.2

# MSnbase 1.7

## Changes in version 1.7.13

 - use donttest instead of dontrun [2015-04-09 Thu]

## Changes in version 1.7.12

 - don't run knntl example to reduce checking time and timeout on
   Windows [2015-04-08 Wed]

## Changes in version 1.7.11

 - fix splitTh, closes issue #49 [2015-04-06 Mon]

## Changes in version 1.7.10


 - Change in vignette to work on zin1: using 12 random th rows. See
 issue #49 for details. Fixed in version 1.7.11. [2015-04-03 Fri]

## Changes in version 1.7.9

 - updated tl vignette [2015-04-02 Thu]
 - depending on latest (1.5.8) pRolocdata [2015-04-02 Thu]

## Changes in version 1.7.8

 - updating tl vig [2015-04-02 Thu]
 - update getParams documentation [2015-04-02 Thu]

## Changes in version 1.7.7

 - renaming theta.scores to knntl.scores [2015-03-24 Tue]

## Changes in version 1.7.6

 - added the theta inductive transfer infrastructure [2015-02-05 Thu]
 - theta vignette stub [2015-02-06 Fri]
 - rename getClasses to getMarkerClasses [2015-02-06 Fri]
 - added the infrastructure to create GO MSnSet [2015-02-07 Sat]
 - Fixed ml vignette [2015-02-16 Mon]
 - new filterZeroRows function [2015-03-10 Tue]
 - hpa data section [2015-03-10 Tue]
 - theta sections [2015-03-10 Tue]
 - deprecate getRegulari[z|s]edParams [2015-03-11 Wed]

## Changes in version 1.7.5

 - Fix vignettes: run bibtex and pdflatex twice and typo
   [2015-02-03 Tue]

## Changes in version 1.7.4

 - Use default Sweave call to build vignette [2015-01-24 Sat]

## Changes in version 1.7.3

 - use Biocstyle [2015-01-23 Fri]
 - replace library/require by requireNamespace [2015-01-23 Fri]

## Changes in version 1.7.2

 - added t-SNE method to plot2D [2015-01-14 Wed]
 - Updated NAMESPACE imports [2015-01-14 Wed]

## Changes in version 1.7.1

 - updated vignettes with markers.orig [2014-10-30 Thu]
 - updated ml tests [2014-10-30 Thu]

## Changes in version 1.7.0

 - new devel version, Bioc 3.1

# MSnbase 1.6

## Changes in version 1.6.0

 - new stable version, Bioc 3.0

# MSnbase 1.5

## Changes in version 1.5.19

 - added Video tag in DESCRIPTION [2014-10-07 Tue]

## Changes in version 1.5.18

 - HUPO 2014 poster [2014-01-02 Thu]

## Changes in version 1.5.17

 - fix 'replacing previous import by MLInterfaces::plot when loading
   pRoloc' warning by using specific importFrom [2014-09-27 Sat]

## Changes in version 1.5.16

 - new pRolocGUI section [2014-08-15 Fri]
 - new foi section [2014-08-16 Sat]

## Changes in version 1.5.15

 - svmOpt sigma defaults changed from 10^(-2:3) to 10^(-3:2)
   [2014-08-15 Fri]
 - in xxxOptimisation, the best parameter(s) for the validation
   classification runs are now chosen at random instead of using the
   first best param (see change in pRoloc:::getBestParam that got a
   sample argument defaulted to TRUE) [2014-08-15 Fri]

 - When calculating macroF1 scores (xval and validation), NAs are set
   to 0 (via MLInterfaces:::.macroF1(..., naAs0 = TRUE)). The macro F1
   will not be NA (when mean of F1s is calculated) but lowered. This
   avoids having an NA macro F1 when 1 (or more) classe(s) end(s) up
   with NA (also set to 0) precision(s) or recall(s) [2014-08-15 Fri]

## Changes in version 1.5.14

 - add title to plotDist figures [2014-08-13 Wed]

## Changes in version 1.5.13

 - none

## Changes in version 1.5.12

 - support mirrorX/mirrorY in highlightOnPlot [2014-07-22 Tue]

## Changes in version 1.5.11

 - Remove plot2D outliers param [2014-07-10 Thu]
 - fix subsetAsDataFrame for keepColNames and write unit test
   [2014-07-11 Fri]

## Changes in version 1.5.10

 - alias lopims4 function [2014-07-01 Tue]

## Changes in version 1.5.9

 - Export single steps of lopims [2014-06-23 Mon]
 - Rephrase classifier parameter optimisation proceudres [2014-06-23 Mon]

## Changes in version 1.5.8

 - nndistx_matrix function added to allow use of query matrix
   when calculating knn distances [2014-06-18 Wed]
 - nndist[x]_[matrix|msnset] are now available using the experted
   nndist method [2014-06-18 Wed]
 - markerSet and unknownSet renamed to markerMSnSet and unknownMSnSet
   [2014-06-19 Thu]
 - functions sampleMSnSet and testMSnSet added [2014-06-19 Thu]
 - fix keepColNames in pRoloc:::subsetAsDataFrame - fcol was
   always renamed to "markers" [2014-06-19 Thu]

## Changes in version 1.5.7

 - add recommended biocView [2014-06-05 Thu]

## Changes in version 1.5.6

 - addMarkers has a new mcol argument to set the markers feature
   variable label [2014-05-29 Thu]

## Changes in version 1.5.5

 - replaced MSVBAR::rmultnorm with mvtnorm::rmvnorm since the former
   has been removed from CRAN and don't import [2014-05-21 Wed]
 - Bug tracking [2014-05-26 Mon]

## Changes in version 1.5.4

 - testMarkers gets an error argument [2014-05-14 Wed]
 - plotDist now has a ylim argument [2014-05-21 Wed]

## Changes in version 1.5.3

 - import all MLInterfaces [2014-04-30 Wed]
 - new param optim secion in ml vignette [2014-05-05 Mon]
 - various ml typos and pRolocmakers man update [2014-05-06 Tue]

## Changes in version 1.5.2

 - In plotDist, ... is now passed to matlines instead of plot and has
   a new lty parameter [2014-04-17 Thu]
 - new highlightOnPlot function, using the new features of interest
   infrastructure [2014-04-29 Tue]

## Changes in version 1.5.1

 - new dunkley2006 pdunit object created with mclust 4.3
   [2014-04-08 Tue]
 - addMarker also accepts fcol and addMarkers unit test
   [2014-04-14 Mon]

## Changes in version 1.5.0

 - Bioc devel 3.0

# MSnbase 1.4

## Changes in version 1.4.0

 - Bioc release 2.14

# MSnbase 1.3

## Changes in version 1.3.19

 - fixed error introduced with mclust 4.3 (that now returns the data
   in the Mclust output - see comment in pRoloc:::gmmOutliers for
   details) [2014-04-07 Mon]

## Changes in version 1.3.18

 - getPredictions can take class-specific scores [2014-04-04 Fri]

## Changes in version 1.3.17

 - fixed newly introduced bug (see 1.3.16) in
   pRoloc:::subsetAsDataFrame - thank you unit tests for saving me,
   again [2014-03-26 Wed]

## Changes in version 1.3.16

 - pRoloc:::subsetAsDataFrame now preserved original sample/column
   names [2014-03-24 Mon]
 - fixed wrong message when using col and pch in plot2D
   [2014-03-25 Tue]

## Changes in version 1.3.15

 - updated pRolocmarkers("mmus") [2014-03-21 Fri]
 - moved extdata/*csv to pRolocdata [2014-03-23 Sun]
 - using *csv from pRolocdata [2014-03-23 Sun]

## Changes in version 1.3.14

 - deleted tab character [2014-03-15 Sat]
 - added support for GMM parametrisation to phenoDisco [2014-03-17 Mon]
 - message instead of warning when using colour and pch [2014-03-21 Fri]
 - remove 1 duplicated mouse marker [2014-03-21 Fri]

## Changes in version 1.3.13

 - fixed a bug in addLegend [2014-03-14 Fri]
 - updated testing to testthat 0.8 [2014-03-14 Fri]
 - Fixing several warnings about symnbols being replaced upon pRoloc
   loading and note about usage of ::: [2014-03-15 Sat]

## Changes in version 1.3.12

 - added phenoDisco2 for testing, allows choice of GMM
   parameters [2014-02-26 Wed]
 - removed duplicated fly markers [2014-02-28 Fri]
 - updated affiliations in vignettes [2014-03-10 Mon]

## Changes in version 1.3.11

 - modify to new biocViews to DESCRIPTION file by s.arora [2014-03-04]

## Changes in version 1.3.10

 - new phenoDisco ndims argument to use more than two two principal
   components as input for discovery analysis [2014-01-03 Mon]
 - fixed and updated phenoDisco logging [2014-02-03 Mon]
 - added support for parallel phenoDisco execution. See BPPARAM
   argument [2014-02-03 Mon]
 - fixed issues with using PD and ndims [2014-02-10 Mon]
 - fixed call to anyUnknown in PD code [2014-02-10 Mon]
 - checking if duplicated markers in addMarkers [2014-02-14 Fri]

## Changes in version 1.3.9

 - bump to force rebuild for new Rcpp

## Changes in version 1.3.8

 - Removed trailing space in mmus nucleus markers [2014-01-21 Tue]
 - using filterNA to remove features with missing values
   in plot2D [2014-01-21 Tue]
 - fixed plot2D/addLegend [2014-01-23 Thu]

## Changes in version 1.3.7

 - fixed addLegend to use correct colours (order) [2014-01-20 Mon]
 - fix typo in addMarkers man [2014-01-20 Mon]
 - re-arranged stockpch so that interleave full and empty
   plotting character [2014-01-20 Mon]
 - Removed last stockcol (tomato), too cose to "#FF7F00" [2014-01-20 Mon]

## Changes in version 1.3.6

 - updated human markers: keep new pd.markers phenotypes,
   validated by Lisa and remove singletons [2014-01-14 Tue]
 - first stockpch is noe 19 [2014-01-16 Thu]
 - plot2D(.., pch) now taken into account for labelled
   data [2014-01-16 Thu]
 - removed alpha plot2D argument [2014-01-17 Fri]
 - updated plot2D and addLegend function with support for more
   organelle groups than colours. The previous versions are
   available as plot2D_v1 and addLegend_v1. [2014-01-17 Fri]

## Changes in version 1.3.5

 - pRoloc citation [2014-01-12 Sun]

## Changes in version 1.3.4

 - new unknownSet function [2013-12-11 Wed]
 - updated vignettes to account for tan2009r1 changes [2013-12-16 Mon]

## Changes in version 1.3.3

 - update citation in phenoDisco.Rd [2013-11-22 Fri]
 - typos in vignettes [2013-11-25 Mon]

## Changes in version 1.3.2

 - new markers [2013-11-06 Wed]
 - markers in vinette [2013-11-15 Fri]
 - mouse markers [2013-11-18 Mon]

## Changes in version 1.3.1

 - using combineFeatures(..., na.rm=TRUE) in lopims pipeine [2013-10-22 Tue]
 - corrected spelling errors in phenoDisco doc [2013-10-29 Tue]

## Changes in version 1.3.0

 - next devel version fot Bioc 2.14

# MSnbase 1.1

## Changes in version 1.1.8

 - semi-sup and sup comparison sections in ml vig [2013-10-09 Wed]

## Changes in version 1.1.7

 - getMarkers has names arg [2013-10-01 Tue]

## Changes in version 1.1.6

 - new outliers arg to plot2D [2013-09-23 Mon]
 - cite addMarkers in vignette [2013-09-26 Thu]
 - add code chunk in poster vig [2013-09-26 Thu]

## Changes in version 1.1.5

 - added biocViews [2013-09-19 Thu]
 - added knitr vig engine in ml [2013-09-19 Thu]
 - import dependencies [2013-09-20 Fri]

## Changes in version 1.1.4

 - building ml vignette in Makefile [2013-09-11 Wed]

## Changes in version 1.1.3

 - new private anyUnknown function, used in phenoDisco, to check
   for the presence of unlabelled ('unknown') data [2013-08-27 Tue]
 - added a note in vignette about "unknown" convention to define
   protein with unknown localisation [2013-08-28 Wed]
 - Added HUPO 2011 poster [2013-09-09 Mon]

## Changes in version 1.1.2

 - fixed Author[s]@R [2013-05-16 Thu]
 - na.rm=TRUE in f1Count [2013-05-19 Sun]
 - added CITATION file [2013-06-07 Fri]
 - new testMarkers function [2013-06-29 Sat]
 - error message in getBestParams suggests to use testMarkers [2013-06-29 Sat]
 - nndist methods (see issue #23) [2013-07-01 Mon]
 - remove unnecessary as.matrix in plot2D's cmdscale [2013-07-19 Fri]
 - added plot2D(..., method = "MDS") example [2013-07-19 Fri]
 - changed 'euclidian' to 'euclidean' in nndist_matrix [2013-07-26 Fri]
 - fixed row ordering in phenoDisco, input and output rownames are now the same [2013-08-13 Tue]
 - Using filterNA in phenoDisco [2013-08-13 Tue]

## Changes in version 1.1.1

 - Update README.md to reflect availability in stable release
   and biocLite installation [2013-04-07 Sun]
 - new MSe pipeline [2013-04-07 Sun]
 - perTurbo using Rcpp [2013-04-16 Tue]
 - initial clustering infrastructure (not exported) [2013-04-17 Wed]
 - new markerSet function [2013-04-24 Wed]
 - plsaOptim's ncomp is now 2:6 [2013-04-27 Sat]
 - added forword to vignette [2013-04-29 Mon]
 - default k is now seq(3, 15, 2) in knnOptim [2013-05-09 Thu]

## Changes in version 1.1.0

 - Bioc 2.13 devel version bump

# MSnbase 1.0

## Changes in version 1.0.0

 - Bioc 2.12 stable version bump

# MSnbase 0.99

## Changes in version 0.99.17

 - illustrating class.weights in the
   vignette [2013-03-24 Sun]

## Changes in version 0.99.16

 - new addMarkers function [2013-03-22 Fri]

## Changes in version 0.99.15

 - Fixing issues in vignette [2013-03-22 Fri]

## Changes in version 0.99.14

 - Added scale in tutorial [2013-03-02 Sat]
 - implemented viction [2013-03-09 Sat]
 - New vignette section on phenoDisco follow
   up classification [2013-03-19 Tue]
 - Using knitr engine [2013-03-19 Tue]

## Changes in version 0.99.13

 - depending on MSnbase ]= 1.7.23 as makeNaData
   needs MSnbase:::nologging [2013-02-27 Wed]
 - updated phenoDisco parameters, added error messages when an
   insufficient number of markers per class and/or number of
   classes are used and updated phenoDisco help file [2012-02-28 Thu]
 - added first belief diffscores function [2013-03-01 Fri]

## Changes in version 0.99.12

 - plot2D has gained a plot argument [2013-02-19 Tue]
 - new f1Count method [2013-02-20 Wed]
 - new makeNaData function [2013-02-20 Wed]
 - new makeNaData2 function [2013-02-21 Thu]
 - new whichNAfunction [2013-02-20 Wed]
 - updated tutorial vignette [2013-02-20 Wed]
 - Now passing ... to predictor functions in
   xxxClassification (reported by Marianne Sandin) [2013-02-26 Tue]

## Changes in version 0.99.11

 - setUnknowncol(NULL) and friedns reset to default
   values [2013-02-19 Tue]

## Changes in version 0.99.10

 - new default col/pch setters [2013-02-16 Sat]

## Changes in version 0.99.9

 - Updates to phenoDisco: verbose param, fixed error in example,
   adding params to processingInfo [2013-02-11 Mon]
 - Unexported getOtherParams method [2013-02-12 Tue]
 - adding export param documentation [2013-02-12 Tue]

## Changes in version 0.99.8

 - Integration of the perTurbo algorithms, contributed
   by Thomas Burger and Samuel Wieczorek [2013-01-18 Fri]
 - summariseMatList now has na.rm = TRUE by default [2013-01-19 Sat]
 - PerTurbo's inv/reg now as other hyperparams [2013-02-11 Mon]

## Changes in version 0.99.7

 - knitr 1.0 compatibility [2013-01-15 Tue]
 - Updated phenoDisco documentation and README [2013-01-10 Thu]

## Changes in version 0.99.6

 - removed updateClass man [2013-01-03 Thu]
 - removed old Rd files [2013-01-03 Thu]
 - Deprecating *Regularisation and *Prediction function [2013-01-03 Thu]
 - Updated new names in test_ml.R [2013-01-04 Fri]

## Changes in version 0.99.5

 - more reg data into GenRegRes objects
   - cmMatrices (knn) [2012-11-15 Thu] (other) [2012-11-30 Fri]
   - testPartitions (knn) [2012-11-17 Sat] (other) [2012-11-30 Fri]
 - new minMarkers function [2012-11-16 Fri]
 - renamed updateClass to minClassScore [2012-11-16 Fri]
 - renamed xxxRegularisation to xxxOptimisation [2012-11-30 Fri]
 - renamed xxxPrediction to xxxClassification [2012-11-30 Fri]
 - renamed getRegularisedParams to getParams [2012-11-30 Fri]
 - moved exprsToRatios to MSnbase [2012-12-05 Wed]
 - updated phenoDisco help file [2012-12-07 Fri]
 - updated phenoDisco code to cope with identical protein profiles [2012-12-07 Fri]

## Changes in version 0.99.4

 - Adding README file describing Rd generation
   and suggesting roxygen2 [2012-11-14 Wed]
 - Added scol=NULL to ignore completely [2012-11-14 Wed]

## Changes in version 0.99.3

 - pdres in extdata - updated vignette [2012-11-14 Wed]
 - udpated pd's GS/times default [2012-11-14 Wed]

## Changes in version 0.99.2

 - fixed MLearn("formula", "MSnSet", "clusteringSchema",
   "missing") - interface was wrong [2012-11-10 Sat]
 - vignette updates [2012-11-10 Sat] [2012-11-11 Sun]
 - nicer knn score names when scores = "all" [2012-11-11 Sun]

## Changes in version 0.99.1

 - updated exprsToRatio when ncol(object) is 2 [2012-11-05 Mon]
 - typos in vignette [2012-11-09 Fri] [2012-11-10 Sat]
 - new MLearn method for signature
   c("formula", "MSnSet", "clusteringSchema", "missing") [2012-11-09 Fri]
 - Several vignette udpates [2012-11-10 Sat]

## Changes in version 0.99.0

 - Submission to Bioc [2012-11-04 Sun]
