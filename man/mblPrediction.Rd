\name{mblPrediction}
\alias{mblPrediction}
\title{mbl prediction}
\usage{
  mblPrediction(object, fcol = "markers", verbose = TRUE,
    scores = c("prediction", "all", "none"), seed, ...)
}
\arguments{
  \item{object}{An instance of class
  code{"\linkS4class{MSnSet}"}.}

  \item{fcol}{The feature meta-data containing marker
  definitions. Default is \code{markers}.}

  \item{verbose}{A \code{logical} defining whether a
  progress bar is displayed.}

  \item{scores}{One of \code{"prediction"}, \code{"all"} or
  \code{"none"}.}

  \item{seed}{The optional random number generator seed.}

  \item{...}{Additional parameters passed to
  \code{\link{svmRegularisation}}.}
}
\value{
  An updated instance of class
  code{"\linkS4class{MSnSet}"}.
}
\description{
  Multiple binary learning algorithm.
}
\author{
  Laurent Gatto
}

