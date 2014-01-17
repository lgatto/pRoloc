\name{addLegend_v1}
\alias{addLegend_v1}
\title{Adds a legend}
\usage{
addLegend_v1(object, fcol = "markers", where = "other", col, ...)
}
\arguments{
  \item{object}{An instance of class \code{MSnSet}}

  \item{fcol}{Feature meta-data label (fData column name)
  defining the groups to be differentiated using different
  colours. Default is \code{markers}.}

  \item{where}{One of \code{"other"}, \code{"bottomleft"},
  \code{"bottomright"}, \code{"topleft"} or
  \code{"topright"} defining the location of the legend.
  \code{"other"} opens a new graphics device, while the
  other locations are passed to \code{\link{legend}}.}

  \item{col}{A \code{character} defining point colours.}

  \item{...}{Additional parameters passed to
  \code{\link{legend}}.}
}
\value{
Invisibly returns \code{NULL}
}
\description{
Adds a legend to a \code{\link{plot2D_v1}} figure.
}
\author{
Laurent Gatto
}

