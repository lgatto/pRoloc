% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/machinelearning-functions-PerTurbo.R
\name{perTurboClassification}
\alias{perTurboClassification}
\title{perTurbo classification}
\usage{
perTurboClassification(
  object,
  assessRes,
  scores = c("prediction", "all", "none"),
  pRegul,
  sigma,
  inv,
  reg,
  fcol = "markers"
)
}
\arguments{
\item{object}{An instance of class \code{"\linkS4class{MSnSet}"}.}

\item{assessRes}{An instance of class
\code{"\linkS4class{GenRegRes}"}, as generated by
\code{\link{svmRegularisation}}.}

\item{scores}{One of \code{"prediction"}, \code{"all"} or
\code{"none"} to report the score for the predicted class
only, for all classes or none.}

\item{pRegul}{If \code{assessRes} is missing, a \code{pRegul} must
be provided.  See \code{\link{perTurboOptimisation}} for
details.}

\item{sigma}{If \code{assessRes} is missing, a \code{sigma} must
be provided.  See \code{\link{perTurboOptimisation}} for
details.}

\item{inv}{The type of algorithm used to invert the matrix.
Values are : "Inversion Cholesky" (\code{\link{chol2inv}}),
"Moore Penrose" (\code{\link{ginv}}), "solve"
(\code{\link{solve}}), "svd" (\code{\link{svd}}).  Default
value is \code{"Inversion Cholesky"}.}

\item{reg}{The type of regularisation of matrix.  Values are
"none", "trunc" or "tikhonov".  Default value is
\code{"tikhonov"}.}

\item{fcol}{The feature meta-data containing marker definitions.
Default is \code{markers}.}
}
\value{
An instance of class \code{"\linkS4class{MSnSet}"} with
    \code{perTurbo} and \code{perTurbo.scores} feature variables
    storing the classification results and scores respectively.
}
\description{
Classification using the PerTurbo algorithm.
}
\examples{
library(pRolocdata)
data(dunkley2006)
## reducing parameter search space 
params <- perTurboOptimisation(dunkley2006,
                               pRegul = 2^seq(-2,2,2),
                               sigma = 10^seq(-1, 1, 1),
                               inv = "Inversion Cholesky",
                               reg ="tikhonov",
                               times = 3)
params
plot(params)
f1Count(params)
levelPlot(params)
getParams(params)
res <- perTurboClassification(dunkley2006, params)
getPredictions(res, fcol = "perTurbo")
getPredictions(res, fcol = "perTurbo", t = 0.75)
plot2D(res, fcol = "perTurbo")
}
\references{
N. Courty, T. Burger,
    J. Laurent. "PerTurbo: a new classification algorithm based on the spectrum perturbations of the Laplace-Beltrami operator",
    The European Conference on Machine Learning and Principles and
    Practice of Knowledge Discovery in Databases (ECML-PKDD 2011),
    D. Gunopulos et al.  (Eds.): ECML PKDD 2011, Part I, LNAI
    6911, pp. 359 - 374, Athens, Greece, September 2011.
}
\author{
Thomas Burger and Samuel Wieczorek
}
