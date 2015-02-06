\name{GenRegRes-class}
\Rdversion{1.1}
\docType{class}
\alias{GenRegRes}
\alias{ThetaRegRes}
\alias{class:GenRegRes}
\alias{class:ThetaRegRes}
\alias{GenRegRes-class}
\alias{getF1Scores}
\alias{getF1Scores,GenRegRes-method}
\alias{f1Count}
\alias{f1Count,GenRegRes-method}
\alias{getParams}
\alias{getParams,GenRegRes-method}
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

\title{Class \code{"GenRegRes"} and \code{"ThetaRegRes"}}

\description{
  Regularisation framework containers.
}

\section{Objects from the Class}{
  Object of this class are created with the respective regularisation
  function: \code{\link{knnOptimisation}},
  \code{\link{svmOptimisation}}, \code{\link{plsdaOptimisation}},
  \code{\link{thetaOptimisation}}, ... 
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

  Only in \code{ThetaRegRes}:

  \describe{
    \item{\code{predictions}:}{ A \code{list} of predictions for the
      optimisation iterations. }    
    \item{\code{otherWeights}:}{ Alternative best theta weigts: a vector
      per iterations, \code{NULL} if no other best weights were found. }    
  }
  
}

\section{Methods}{
  \describe{
    \item{getF1Scores}{Returns a matrix of F1 scores for the
      optimisation parameters. }
    \item{f1Count}{\code{signature(object = "GenRegRes", t =
	"numeric")} and \code{signature(object = "ThetaRegRes", t =
	"numeric")}: Constructs a table of all possible parameter
      combination and count how many have an F1 scores greater or equal
      than \code{t}. When \code{t} is missing (default), the best F1
      score is used. This method is useful in conjunctin with
      \code{plot}. }
    \item{getParams}{Returns the \emph{best} parameters. It is however
      recommended to inspect the optimisation results. }
    \item{getSeed}{Returns the seed used for the optimisation run. }
    \item{getWarnings}{\code{signature(object = "GenRegRes")}: Returns a
      vector of recorded warnings. }
    \item{levelPlot}{\code{signature(object = "GenRegRes")}: Plots a
      heatmap of of the optimisation results. Only for
      \code{"GenRegRes"} instances. }
    \item{plot}{Plots the optisisation results. }
    \item{show}{Shows the object. }
  }
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\examples{
showClass("GenRegRes")
showClass("ThetaRegRes")
}

\keyword{classes}
