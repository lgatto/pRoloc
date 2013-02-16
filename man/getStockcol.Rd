\name{getStockcol}
\alias{getStockcol}
\alias{getStockpch}
\alias{getUnknowncol}
\alias{getUnknownpch}
\alias{setStockcol}
\alias{setStockpch}
\alias{setUnknowncol}
\alias{setUnknownpch}
\title{Manage default colours and point characters}
\usage{
  getStockcol()

  setStockcol(cols)

  getStockpch()

  setStockpch(pchs)

  getUnknowncol()

  setUnknowncol(col)

  getUnknownpch()

  setUnknownpch(pch)
}
\arguments{
  \item{cols}{A vector of colour \code{characters}.}

  \item{pchs}{A vector of \code{numeric}.}

  \item{col}{A colour \code{character}.}

  \item{pch}{A \code{numeric} vector of length 1.}
}
\value{
  A \code{character} vector.

  Invisibly returns \code{cols}.

  A \code{numeric} vector.

  Invisibly returns \code{pchs}.

  A \code{character} vector or length 1.

  Invisibly returns \code{col}.

  A \code{numeric} vector of length 1.

  Invisibly returns \code{pch}.
}
\description{
  These functions allow to get/set the default colours and
  point character that are used when plotting organelle
  clusters and unknown features. These values are
  parametrised at the session level.
}
\examples{
## defaults for clusters
getStockcol()
getStockpch()
## unknown features
getUnknownpch()
getUnknowncol()
## an example
library(pRolocdata)
data(dunkley2006)
par(mfrow = c(2, 1))
plot2D(dunkley2006, fcol = "markers", main = 'Default colours')
setUnknowncol("black")
plot2D(dunkley2006, fcol = "markers", main = 'setUnknowncol("black")')
}
\author{
  Laurent Gatto
}

