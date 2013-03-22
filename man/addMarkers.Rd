\name{addMarkers}
\alias{addMarkers}
\title{Adds markers to the data}
\usage{
  addMarkers(object, markerfile, verbose = TRUE)
}
\arguments{
  \item{object}{An instance of class \code{MSnSet}.}

  \item{markerfile}{A \code{character} with the name the
  markers' csv file.}

  \item{verbose}{A \code{logical} indicating if number of
  markers and marker table should be printed to the
  console.}
}
\value{
  A new instance of class \code{MSnSet} with an additional
  \code{markers} feature variable.
}
\description{
  The function adds a 'markers' feature variable. These
  markers are read from a comma separated values (csv)
  spreadsheet file. This markers file is expected to have 2
  columns (others are ignored) where the first is the name
  of the marker features and the second the group label. It
  is essential to assure that \code{featureNames(object)}
  and marker names (first column) match, i.e. the same
  feature identifiers and case fold are used.
}
\author{
  Laurent Gatto
}

