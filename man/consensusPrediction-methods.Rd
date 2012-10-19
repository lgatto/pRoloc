\name{consensusPrediction-methods}
\docType{methods}
\alias{consensusPrediction}
\alias{consensusPrediction-methods}
\alias{consensusPrediction,MSnSet-method}

\title{ Ensembl method to combine classification results }

\description{
  This method combines a set of classification results, as obtained by
  \code{\link{svmPrediction}}, \code{\link{knnPrediction}},
  \code{\link{nnetPrediction}}, \ldots using one of three methods,
  namely \code{"filteredMajority"}, \code{"majority"} or
  \code{"strict"}. \code{"strict"} requires classification predictions
  to be exactly equal across prediciton methods. \code{"majority"} uses
  a majority voting mechanism, while \code{"filteredMajority"} first
  sets all predictions whose score are lower than the result of
  \code{fun} on that column score to \code{"unknown"}, prior to apply
  majority voting.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnSet", method = "character", fcols
	= "character", fun = "function")}}{
      See description.
    }
  }
}
\keyword{methods}

