\name{getAttributesOfInterest}
\alias{getAttributesOfInterest}
\title{Get Biomart attributes of interest.}
\usage{
  getAttributesOfInterest()
}
\value{
  A \code{character} with the attributes that the biomart
  data sets used in pRoloc all have.
}
\description{
  This function gives the attributes that are used as a
  criterion for ENSEMBL data sets to be used in pRoloc.
  These are the attributes that are susceptible to be
  queried by pRoloc functions as \code{\link{getGO}} or
  \code{\link{getSeq}}.
}
\author{
  Laurent Gatto
}

