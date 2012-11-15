\name{getPredictions}
\alias{getPredictions}
\title{Returns the predictions in an 'MSnSet'}
\usage{
  getPredictions(object, fcol, scol, t = 0, verbose = TRUE)
}
\arguments{
  \item{object}{An instance of class
  \code{"\linkS4class{MSnSet}"}.}

  \item{fcol}{The name of the prediction column in the
  \code{featureData} slot.}

  \item{scol}{The name of the prediction score column in
  the \code{featureData} slot. If missing, created by
  pasting '.scores' after \code{fcol}. If \code{NULL},
  ignored.}

  \item{t}{The score threshold. Predictions with score < t
  are set to 'unknown'. Default is 0.}

  \item{verbose}{If \code{TRUE}, a prediction table is
  printed and the predictions are returned invisibly. If
  \code{FALSE}, the predictions are returned.}
}
\value{
  A \code{character} of length \code{ncol(object)}.
}
\description{
  Convenience accessor to the predicted feature
  localisation in an 'MSnSet'. This function returns the
  predictions of an \code{MSnSet} instance. As a side
  effect, it prints out a prediction table.
}
\author{
  Laurent Gatto
}

