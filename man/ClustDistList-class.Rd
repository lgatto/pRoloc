\name{ClustDistList-class}
\Rdversion{1.1}
\docType{class}
\alias{ClustDistList}
\alias{class:ClustDistList}
\alias{ClustDistList-class}
\alias{plot,ClustDistList,missing-method}
\alias{show,ClustDistList-method}
\alias{[,ClustDistList,ANY,ANY,ANY-method}
\alias{[,ClustDistList,ANY,missing,missing-method}
\alias{[[,ClustDistList,ANY,ANY-method}
\alias{[[,ClustDistList,ANY,missing-method}
\alias{length,ClustDistList-method}
\alias{names,ClustDistList-method}
\alias{names<-,ClustDistList,ANY-method}
\alias{lapply,ClustDistList-method}
\alias{sapply,ClustDistList-method}


\title{Storing multiple ClustDist instances}

\description{
  A class for storing lists of \code{\linkS4class{ClustDist}}
  instances. 
}

\section{Objects from the Class}{
  Object of this class are created with the \code{clustDist}
  function.
}

\section{Slots}{
  \describe{
    
    \item{\code{x}:}{Object of class \code{list} containing valid
      \code{ClustDist} instances.}

    \item{\code{log}:}{Object of class \code{list} containing an object
      creation log, containing among other elements the call that
      generated the object.}

    \item{\code{.__classVersion__}:}{The version of the instance. For
      development purposes only. }
  }
}
  

\section{Methods}{
  \describe{
    \item{\code{"[["}}{Extracts a single \code{ClustDist} at position. }
    \item{\code{"["}}{Extracts one of more \code{ClustDists} as
      \code{ClustDistList}. }
    \item{\code{length}}{Returns the number of \code{ClustDists}. }
    \item{\code{names}}{Returns the names of \code{ClustDists}, if
      available. The replacement method is also available. }
    \item{\code{show}}{ Display the object by printing a short
      summary. }
    \item{\code{lapply(x, FUN, ...)}}{ Apply function \code{FUN} to each
      element of the input \code{x}. If the application of \code{FUN}
      returns and \code{ClustDist}, then the return value is an
      \code{ClustDistList}, otherwise a \code{list}}.
    \item{\code{plot}}{Plots a boxplot of the distance results per
      protein set.}
  }
}


\author{
  Lisa M Breckels <lms79@cam.ac.uk>
}


\examples{
  library('pRolocdata')
  data(dunkley2006)
  par <- setAnnotationParams(inputs =
                    c("Arabidopsis thaliana genes",
                    "Gene stable ID"))
                    
  ## add protein set/annotation information                  
  xx <- addGoAnnotations(dunkley2006, par)
  
  ## filter
  xx <- filterMinMarkers(xx, n = 50)
  xx <- filterMaxMarkers(xx, p = .25)
  
  ## get distances for protein sets
  dd <- clustDist(xx)
  
  ## plot distances for all protein sets 
  plot(dd)
  
  names(dd)
  
  ## Extract first 4 ClustDist objects of the ClustDistList
  dd[1:4]
  
  ## Extract 1st ClustDist object
  dd[[1]]
}

\keyword{classes}
