\name{plot2D}
\alias{plot2D}
\title{Plot organelle assignment data and results.}
\usage{
  plot2D(object, fcol = "markers", fpch,
    unknown = "unknown", dims = 1:2, alpha, score = 1,
    method = c("PCA", "MDS"), axsSwitch = FALSE,
    mirrorX = FALSE, mirrorY = FALSE, col, pch, cex,
    identify = FALSE, ...)
}
\arguments{
  \item{object}{An instance of class \code{MSnSet}.}

  \item{fcol}{Feature meta-data label (fData column name)
  defining the groups to be differentiated using different
  colours. Default is \code{markers}. Use \code{NULL} to
  suppress any colouring.}

  \item{fpch}{Featre meta-data label (fData column name)
  desining the groups to be differentiated using different
  point symbols.}

  \item{unknown}{A \code{character} (default is
  \code{"unknown"}) defining how proteins of unknown
  localisation are labelled.}

  \item{dims}{A \code{numeric} of length 2 defining the
  dimensions to be plotted, i.e the PC/MDS axes.}

  \item{alpha}{A numeric defining the alpha channel
  (transparency) of the points, where \code{0 <= alpha <=
  1}, 0 and 1 being completely transparent and opaque.}

  \item{score}{A numeric specifying the minimum organelle
  assignment score to consider features to be assigned an
  organelle. (not yet implemented).}

  \item{method}{One of \code{PCA} (default) or \code{MDS},
  defining if dimensionality reduction is done using
  principal component analysis (see \code{\link{prcomp}})
  or classical multidimensional scaling (see
  \code{\link{cmdscale}}).}

  \item{axsSwitch}{A \code{logical} indicating whether the
  axes should be switched.}

  \item{mirrorX}{A \code{logical} indicating whether the x
  axis should be mirrored?}

  \item{mirrorY}{A \code{logical} indicating whether the y
  axis should be mirrored?}

  \item{col}{A \code{character} of appropriate length
  defining colours.}

  \item{pch}{A \code{character} of appropriate length
  defining point character.}

  \item{cex}{Character expansion.}

  \item{identify}{A logical (default is \code{TRUE})
  defining if user interaction will be expected to identify
  individual data points on the plot. See also
  \code{\link{identify}}.}

  \item{...}{Additional parameters passed to \code{plot}
  and \code{points}.}
}
\value{
  Used for its side effects of generating a plot. Invisibly
  returns the 2d data.
}
\description{
  Generate 2 dimensional or feature distribution plots to
  illustrate localistation clusters. In \code{plot2D}, rows
  containing \code{NA} values are removed prior to
  dimention reduction.
}
\examples{
library("pRolocdata")
data(dunkley2006)
plot2D(dunkley2006, fcol = NULL)
plot2D(dunkley2006, fcol = "markers")
addLegend(dunkley2006,
          fcol = "markers",
          where = "topright",
          cex = 0.5, bty = "n", ncol = 3)
title(main = "plot2D example")
}
\author{
  Laurent Gatto <lg390@cam.ac.uk>
}
\seealso{
  \code{\link{addLegend}} to add a legend to \code{plot2D}
  figures and \code{\link{plotDist}} for alternative
  graphical representation of quantitative organelle
  proteomics data.
}

