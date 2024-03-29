% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/machinelearning-functions-nnet.R
\name{nnetClassification}
\alias{nnetClassification}
\alias{nnetPrediction}
\title{nnet classification}
\usage{
nnetClassification(
  object,
  assessRes,
  scores = c("prediction", "all", "none"),
  decay,
  size,
  fcol = "markers",
  ...
)
}
\arguments{
\item{object}{An instance of class \code{"\linkS4class{MSnSet}"}.}

\item{assessRes}{An instance of class
\code{"\linkS4class{GenRegRes}"}, as generated by
\code{\link{nnetOptimisation}}.}

\item{scores}{One of \code{"prediction"}, \code{"all"} or
\code{"none"} to report the score for the predicted class
only, for all classes or none.}

\item{decay}{If \code{assessRes} is missing, a \code{decay} must
be provided.}

\item{size}{If \code{assessRes} is missing, a \code{size} must be
provided.}

\item{fcol}{The feature meta-data containing marker definitions.
Default is \code{markers}.}

\item{...}{Additional parameters passed to \code{\link{nnet}} from
package \code{nnet}.}
}
\value{
An instance of class \code{"\linkS4class{MSnSet}"} with
    \code{nnet} and \code{nnet.scores} feature variables storing
    the classification results and scores respectively.
}
\description{
Classification using the artificial neural network
algorithm.
}
\examples{
library(pRolocdata)
data(dunkley2006)
## reducing parameter search space and iterations 
params <- nnetOptimisation(dunkley2006, decay = 10^(c(-1, -5)), size = c(5, 10), times = 3)
params
plot(params)
f1Count(params)
levelPlot(params)
getParams(params)
res <- nnetClassification(dunkley2006, params)
getPredictions(res, fcol = "nnet")
getPredictions(res, fcol = "nnet", t = 0.75)
plot2D(res, fcol = "nnet")
}
\author{
Laurent Gatto
}
