\name{knnOptimisation}
\alias{knnOptimisation}
\alias{knnOptimization}
\alias{knnRegularisation}
\title{knn parameter optimisation}
\usage{
  knnOptimisation(object, fcol = "markers", k = 3:12,
    times = 100, test.size = 0.2, xval = 5, fun = mean,
    seed, verbose = TRUE, ...)
}
\arguments{
  \item{object}{An instance of class
  \code{"\linkS4class{MSnSet}"}.}

  \item{fcol}{The feature meta-data containing marker
  definitions. Default is \code{markers}.}

  \item{k}{The hyper-parameter. Default values are
  \code{3:12}.}

  \item{times}{The number of times internal
  cross-validation is performed. Default is 100.}

  \item{test.size}{The size of test data. Default is 0.2
  (20 percent).}

  \item{xval}{The \code{n}-cross validation. Default is 5.}

  \item{fun}{The function used to summarise the \code{xval}
  macro F1 matrices.}

  \item{seed}{The optional random number generator seed.}

  \item{verbose}{A \code{logical} defining whether a
  progress bar is displayed.}

  \item{...}{Additional parameters passed to
  \code{\link{knn}} from package \code{class}.}
}
\value{
  An instance of class \code{"\linkS4class{GenRegRes}"}.
}
\description{
  Classification parameter optimisation for the k-nearest
  neighbours algorithm.
}
\author{
  Laurent Gatto
}
\seealso{
  \code{\link{knnClassification}} and example therein.
}

