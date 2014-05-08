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
This function retrieves a list of organelle markers or, if
no \code{species} is provided, prints a description of
available marker sets. The markers can be added to and
\code{MSnSet} using the \code{\link{addMarkers}} function.
}
\details{
The markers have been contributed by various members of the
Cambridge Centre for Proteomics, in particular Dan
Nightingale for yeast, Dr Andy Christoforou for human, Dr
Arnoud Groen for Arabodopsis and Dr Claire Mulvey for
mouse.  In addition, original (curated) markers from the
\code{pRolocdata} datasets have been extracted (see
\code{pRolocdata} for details and references). Curation
involved verification of publicly available subcellular
localisation annotation based on the curators knowledge of
the organelles/proteins considered and tracing the original
statement in the literature.

These markers are provided as a starting point to generate
reliable sets of organelle markers but still need to be
verified against any new data in the light of the
quantitative data and the study conditions.
}
\examples{
pRolocmarkers()
table(pRolocmarkers("atha"))
table(pRolocmarkers("hsap"))
}
\author{
Laurent Gatto
}

