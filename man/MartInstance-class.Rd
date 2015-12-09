\name{MartInstance-class}
\Rdversion{1.1}
\docType{class}
\alias{MartInstance-class}
\alias{MartInstance}
\alias{show,MartInstance-method}

\alias{MartInstanceList-class}
\alias{MartInstanceList}

\alias{as.data.frame.MartInstanceList}
\alias{as.data.frame.MartInstance}

\alias{[,MartInstanceList-method}
\alias{[,MartInstanceList,ANY,ANY-method}
\alias{[,MartInstanceList,ANY,ANY,ANY-method}
\alias{[[,MartInstanceList-method}
\alias{[[,MartInstanceList,ANY,ANY-method}

\alias{sapply,MartInstanceList-method}
\alias{sapply,MartInstanceList,ANY-method}
\alias{lapply,MartInstanceList-method}
\alias{lapply,MartInstanceList,ANY-method}

\alias{nDatasets}
\alias{filterAttrs}
\alias{getMartInstanceList}
\alias{getMartTab}
\alias{getFilterList}

\title{Class \code{"MartInstance"}}

\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}

\section{Objects from the Class}{
  Objects can be created with the respective constructors
  \code{MartInstance} and \code{MartInstanceList}.

}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{host}:}{Object of class \code{"character"} ~~ }
    \item{\code{path}:}{Object of class \code{"character"} ~~ }
    \item{\code{datasets}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{created}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "MartInstance")}: ... }
	 }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

