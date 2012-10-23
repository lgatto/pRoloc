\name{rfRegularisation}
\alias{rfRegularisation}
\title{svm regularisation}
\usage{
  rfRegularisation(object, fcol = "markers", mtry = NULL,
    times = 100, test.size = 0.2, xval = 5, fun = mean,
    seed, verbose = TRUE, ...)
}
\arguments{
  \item{object}{An instance of class
  \code{"\linkS4class{MSnSet}"}.}

  \item{fcol}{The feature meta-data containing marker
  definitions. Default is \code{markers}.}

  \item{mtry}{The hyper-parameter. Default value is
  \code{NULL}.}

  \item{times}{The number of times internal
  cross-validation is performed. Default is 100.}

  \item{test.size}{The size of test data. Default is 0.2
  (20 percent).}

  \item{xval}{The \code{n}-cross validation. Default is 5.}

  \item{fun}{The function used to summarise the
  \code{times} macro F1 matrices.}

  \item{seed}{The optional random number generator seed.}

  \item{verbose}{A \code{logical} defining whether a
  progress bar is displayed.}

  \item{...}{Additional parameters passed to
  \code{\link{randomForest}} from package
  \code{randomForest}.}
}
\value{
  An instance of class \code{"\linkS4class{GenRegRes}"}.
}
\description{
  Hyper-parameters regularisation for the random forest
  algorithm.
}
\author{
  Laurent Gatto
}

