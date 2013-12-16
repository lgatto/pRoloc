\name{svmOptimisation}
\alias{svmOptimisation}
\alias{svmOptimization}
\alias{svmRegularisation}
\title{svm parameter optimisation}
\usage{
svmOptimisation(object, fcol = "markers", cost = 2^(-4:4),
  sigma = 10^(-2:3), times = 100, test.size = 0.2, xval = 5,
  fun = mean, seed, verbose = TRUE, ...)
}
\arguments{
  \item{object}{An instance of class
  \code{"\linkS4class{MSnSet}"}.}

  \item{fcol}{The feature meta-data containing marker
  definitions. Default is \code{markers}.}

  \item{cost}{The hyper-parameter. Default values are
  \code{2^-4:4}.}

  \item{sigma}{The hyper-parameter. Default values are
  \code{10^(-2:3)}.}

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
  \code{\link{svm}} from package \code{e1071}.}
}
\value{
An instance of class \code{"\linkS4class{GenRegRes}"}.
}
\description{
Classification parameter optimisation for the support
vector machine algorithm.
}
\author{
Laurent Gatto
}
\seealso{
\code{\link{svmClassification}} and example therein.
}

