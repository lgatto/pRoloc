\name{GenRegRes-class}
\Rdversion{1.1}
\docType{class}
\alias{GenRegRes}
\alias{ThetaRegRes}
\alias{class:GenRegRes}
\alias{class:ThetaRegRes}
\alias{GenRegRes-class}
\alias{ThetaRegRes-class}
\alias{getF1Scores}
\alias{getF1Scores,GenRegRes-method}
\alias{getF1Scores,ThetaRegRes-method}
\alias{f1Count}
\alias{f1Count,GenRegRes-method}
\alias{f1Count,ThetaRegRes-method}
\alias{getParams}
\alias{getParams,GenRegRes-method}
\alias{getParams,ThetaRegRes-method}
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
\alias{plot,ThetaRegRes,missing-method}
\alias{show,GenRegRes-method}
\alias{show,ThetaRegRes-method}

\alias{combineThetaRegRes}
\alias{favourPrimary}

\title{Class \code{"GenRegRes"} and \code{"ThetaRegRes"}}

\description{
  Regularisation framework containers.
}

\section{Objects from the Class}{
  Object of this class are created with the respective regularisation
  function: \code{\link{knnOptimisation}},
  \code{\link{svmOptimisation}}, \code{\link{plsdaOptimisation}},
  \code{\link{knntlOptimisation}}, ... 
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

    \item{getParams}{

      Returns the \emph{best} parameters. It is however strongly
      recommended to inspect the optimisation results. For a
      \code{ThetaRegRes} optimisation result, the method to chose the
      best parameters can be \code{"median"} (default) or \code{"mean"}
      (the median or mean of the best weights is chosen), \code{"max"}
      (the first weights with the highest macro-F1 score, considering
      that multiple max scoring combinations are possible) or
      \code{"count"} (the observed weight that get the maximum number of
      observations, see \code{f1Count}). The \code{favourP} argument can
      be used to prioritise weights that favour the primary data
      (i.e. heigh weights). See \code{favourPrimary} below.

    }
    
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


\section{Other functions}{
 
  Only for \code{ThetaRegRes}:

  \describe{

    \item{\code{combineThetaRegRes(object)}}{Takes a \code{list} of
      \code{ThetaRegRes} instances to be combined and returnes a new
      \code{ThetaRegRes} instance. }

    \item{\code{favourPrimary(primary, auxiliary, object, verbose =
	TRUE)}}{Takes the \code{primary} and \code{auxiliary} data
	sources (two \code{\linkS4class{MSnSet}} instances) and a
	\code{ThetaRegRes} object and returns and updated
	\code{ThetaRegRes} instance containing best parameters/weigths
	(see the \code{getParams} function) favouring the primary data
	when multiple best theta weights are available. }

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
