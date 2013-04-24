\name{markerSet}
\alias{markerSet}
\title{Extract marker subset}
\usage{
  markerSet(object, fcol = "markers")
}
\arguments{
  \item{object}{An instance of class \code{MSnSet}}

  \item{fcol}{The name of the feature data column, that
  will be used to separate the markers from the proteins of
  unknown localisation (with \code{fData(object)[, fcol] ==
  "unknown")}). Default is to use \code{"markers"}.}
}
\value{
  An new \code{MSnSet} with marker proteins only.
}
\description{
  This function extracts the marker proteins into a new
  \code{MSnSet}.
}
\examples{
library("pRolocdata")
data(dunkley)
mrk <- markerSet(dunkley2006)
dim(dunkley2006)
dim(mrk)
table(fData(dunkley2006)$markers)
table(fData(mrk)$markers)
}
\author{
  Laurent Gatto
}

