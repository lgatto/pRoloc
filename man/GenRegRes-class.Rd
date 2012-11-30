\name{GenRegRes-class}
\Rdversion{1.1}
\docType{class}
\alias{GenRegRes}
\alias{class:GenRegRes}
\alias{GenRegRes-class}
\alias{getF1Scores}
\alias{getF1Scores,GenRegRes-method}
\alias{getRegularisedParams}
\alias{getRegularisedParams,GenRegRes-method}
\alias{getRegularizedParams}
\alias{getRegularizedParams,GenRegRes-method}
\alias{getSeed}
\alias{getSeed,GenRegRes-method}
\alias{getWarnings}
\alias{getWarnings,GenRegRes-method}
\alias{levelPlot}
\alias{levelPlot,GenRegRes-method}
\alias{plot,GenRegRes,missing-method}
\alias{show,GenRegRes-method}

\title{Class \code{"GenRegRes"}}

\description{
  Regularisation framework container.
}

\section{Objects from the Class}{
  Object of this class are created with the respective regularisation
  function: \code{\link{knnRegularisation}},
  \code{\link{svmRegularisation}}, \code{\link{plsdaRegularisation}},
  ... 
}

\section{Slots}{
  \describe{
    \item{\code{algorithm}:}{Object of class \code{"character"} storing
      the machine learning algorithm name. }
    \item{\code{hyperparameters}:}{Object of class \code{"list"} with
      the respective algorithm hyper-parameters tested. }
    \item{\code{design}:}{Object of class \code{"numeric"} describing
      the cross-validation design, the test data size and the number of
      replications.  }
    \item{\code{log}:}{Object of class \code{"list"} with warnings
      thrown during the hyper-parameters regularisation. }
    \item{\code{seed}:}{Object of class \code{"integer"} with the random
      number generation seed. }
    \item{\code{results}:}{Object of class \code{"matrix"} of
      dimenstions \code{times} (see \code{design}) by number of
      hyperparameters + 1 storing the macro F1 values for the respective
      best hyper-parameters for each replication. }
    \item{\code{f1Matrices}:}{Object of class \code{"list"} with
      respective \code{times} cross-validation F1 matrices.  }
    \item{\code{cmMatrices}:}{Object of class \code{"list"} with
      respective \code{times} contingency matrices. }
    \item{\code{testPartitions}:}{Object of class \code{"list"} with
      respective \code{times} test partitions. }    
    \item{\code{datasize}:}{Object of class \code{"list"} with details
      about the respective inner and outter training and testing data
      sizes.  }
  }
}
\section{Methods}{
  \describe{
    \item{getF1Scores}{\code{signature(object = "GenRegRes")}: ... }
    \item{getRegularisedParams}{\code{signature(object = "GenRegRes")}: ... }
    \item{getRegularizedParams}{\code{signature(object = "GenRegRes")}: ... }
    \item{getSeed}{\code{signature(object = "GenRegRes")}: ... }
    \item{getWarnings}{\code{signature(object = "GenRegRes")}: ... }
    \item{levelPlot}{\code{signature(object = "GenRegRes")}: ... }
    \item{plot}{\code{signature(x = "GenRegRes", y = "missing")}: ... }
    \item{show}{\code{signature(object = "GenRegRes")}: ... }
  }
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\examples{
showClass("GenRegRes")
}

\keyword{classes}
