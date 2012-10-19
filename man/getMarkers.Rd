\name{getMarkers}
\alias{getMarkers}
\title{Returns the organelle markers in an 'MSnSet'}
\usage{
  getMarkers(object, fcol = "markers", verbose = TRUE)
}
\arguments{
  \item{object}{An instance of class \code{MSnSet}.}

  \item{fcol}{The name of the markers column in the
  \code{featureData} slot. Default is \code{markers}.}

  \item{verbose}{If \code{TRUE}, a marker table is printed
  and the markers are returned invisibly. If \code{FALSE},
  the markers are returned.}
}
\value{
  A \code{character} of length \code{ncol(object)}.
}
\description{
  Convenience accessor to the organelle markers in an
  'MSnSet'. This function returns the organelle markers of
  an \code{MSnSet} instance. As a side effect, it print out
  a marker table.
}
\examples{
library("pRolocdata")
data(dunkley2006)
mymarkers <- getMarkers(dunkley2006)
}
\author{
  Laurent Gatto
}

