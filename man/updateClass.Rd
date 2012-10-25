\name{updateClass}
\alias{updateClass}
\title{Updates classes based on prediction scores}
\usage{
  updateClass(object, fcol, scol, t = 0)
}
\arguments{
  \item{object}{An instance of class
  \code{"\linkS4class{MSnSet}"}.}

  \item{fcol}{The name of the markers column in the
  \code{featureData} slot.}

  \item{scol}{The name of the prediction score column in
  the \code{featureData} slot. If missing, created by
  pasting '.scores' after \code{fcol}.}

  \item{t}{The score threshold. Predictions with score < t
  are set to 'unknown'. Default is 0.}
}
\value{
  The original \code{object} with a modified
  \code{fData(object)[, fcol]} feature variable.
}
\description{
  This functions updates the classification results in an
  \code{"\linkS4class{MSnSet}"} based on a prediction score
  threshold \code{t}. All features with a score < t are set
  to unknown. Note that the original levels are preserved.
}
\examples{
library(pRolocdata)
data(dunkley)
## random scores
fData(dunkley2006)$assigned.scores <- runif(nrow(dunkley2006))
getPredictions(dunkley2006, fcol = "assigned")
getPredictions(dunkley2006, fcol = "assigned", t = 0.5)
x <- updateClass(dunkley2006, fcol = "assignes", t = 0.5)
getPredictions(x, fcol = "assigned")
all.equal(getPredictions(dunkley2006, fcol = "assigned", t = 0.5),
          getPredictions(x, fcol = "assigned"))
}
\author{
  Laurent Gatto
}

