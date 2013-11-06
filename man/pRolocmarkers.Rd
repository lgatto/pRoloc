\name{pRolocmarkers}
\alias{pRolocmarkers}
\title{Organelle markers}
\usage{
  pRolocmarkers(species)
}
\arguments{
  \item{species}{The species of interest.}
}
\value{
  Prints a description of the available marker lists if
  \code{species} is missing or a named character with
  organelle markers.
}
\description{
  This function retrieves a list of organelle markers or,
  if not \code{species} is provided, prints a description
  of available markers sets. The markers can be added to
  and \code{MSnSet} using the \code{\link{addMarkers}}
  function.
}
\examples{
pRolocmarkers()
table(pRolocmarkers("atha"))
table(pRolocmarkers("hsap"))
}
\author{
  Laurent Gatto
}

