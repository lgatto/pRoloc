\name{makeNaData}
\alias{makeNaData}
\alias{whichNA}
\title{Create a data with missing values}
\usage{
  makeNaData(object, nNA, pNA)

  whichNA(x)
}
\arguments{
  \item{object}{An instance of class \code{MSnSet}.}

  \item{nNA}{The absolute number of missing values to be
  assigned.}

  \item{pNA}{The proportion of missing values to be
  assignmed.}

  \item{x}{A \code{matrix} or an instance of class
  \code{MSnSet}.}
}
\value{
  An instance of class \code{MSnSet}, as \code{x}, but with
  the appropriate number/proportion of missing values. The
  returned object has an additional feature meta-data
  columns, \code{nNA}
}
\description{
  This function takes an instance of class
  \code{"\linkS4class{MSnSet}"} and sets randomly selected
  values to NA. The \code{whichNA} can be used to extract
  the indices of the missing values, as illustrated in the
  example.
}
\examples{
library(pRolocdata)
data(dunkley2006)
sum(is.na(dunkley2006))
dunkleyNA <- makeNaData(dunkley2006, nNA = 150)
sum(is.na(dunkleyNA))
naIdx <- whichNA(dunkleyNA)
head(naIdx)
}
\author{
  Laurent Gatto
}

