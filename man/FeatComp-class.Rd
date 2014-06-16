\name{FeatComp-class}
\Rdversion{1.1}
\docType{class}
\alias{FeatComp-class}
\alias{compfnames,FeatComp,FeatComp-method}
\alias{compfnames}
\alias{show,FeatComp-method}

\title{Class \code{"FeatComp"}}

\description{
  Comparing feature names of two comparable \code{MSnSet} instances.
}

\section{Objects from the Class}{

  Objects can be created with \code{compfnames}. The method compares the
  feature names of two objects of class \code{"MSnSet"}. It prints a
  summary matrix of common and unique feature names and invisibly
  returns a list of \code{FeatComp} instances.

  The function will compute the common and unique features for all
  feature names of the two input objects (\code{featureNames(x)} and
  \code{feautreNames(y)}) as well as distinct subsets as defined in the
  \code{fcol1} and \code{fcol2} feautre variables.

}

\section{Slots}{
  \describe{

    \item{\code{name}:}{Object of class \code{"character"} defining the
      name of the compared features. By convention, \code{"all"} is used
      when all feature names are used; otherwise, the respective levels of
      the feature variables \code{fcol1} and \code{fcol2}. }

    \item{\code{common}:}{Object of class \code{"character"} with the common
      feature names. }

    \item{\code{unique1}:}{Object of class \code{"character"} with the
      features unique to the first \code{MSnSet} (\code{x} in
      \code{comfname}). }

    \item{\code{unique1}:}{Object of class \code{"character"} with the
      features unique to the seconn \code{MSnSet} (\code{y} in
      \code{comfname}). }

    \item{\code{all}:}{Object of class \code{"logical"} defining if all
      features of only a subset were compared. One expects that
      \code{name == "all"} when \code{all} is \code{TRUE}. }
    
  }
}

\section{Methods}{
  \describe{

    \item{compfnames}{\code{signature(x = "FeatComp", y = "FeatComp")}:
      creates the \code{FeatComp} comparison object for instances
      \code{x} and \code{y}. }

    \item{show}{\code{signature(object = "FeatComp")}: prints a summary
      of the \code{object}. }

  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk> and Thomas Naake
}


\examples{
library("pRolocdata")
data(tan2009r1)
data(tan2009r2)
x <- compfnames(tan2009r1, tan2009r2)
x[[1]]
x[2:3]
}

\keyword{classes}
