\name{ClustDist-class}
\Rdversion{1.1}
\docType{class}
\alias{ClustDist}
\alias{class:ClustDist}
\alias{ClustDist-class}
\alias{plot,ClustDist,MSnSet-method}
\alias{show,ClustDist-method}

\title{Class \code{"ClustDist"}}

\description{
  The \code{ClustDist} summaries algorithm information, from 
  running the \code{clustDist} function, such as the number 
  of k's tested for the kmeans, and mean and normalised 
  pairwise (Euclidean) distances per numer of component 
  clusters tested. 
}

\section{Objects from the Class}{
  Object of this class are created with the \code{clustDist}
  function.
}

\section{Slots}{
  \describe{
    \item{\code{k}:}{Object of class \code{"numeric"} storing
      the number of k clusters tested. }
    \item{\code{dist}:}{Object of class \code{"list"} storing 
      the list of distance matrices. }
    \item{\code{term}:}{Object of class \code{"character"} describing
      GO term name.}
    \item{\code{id}:}{Object of class \code{"character"} describing
      the GO term ID.}
    \item{\code{nrow}:}{Object of class \code{"numeric"} showing
      the number of instances in the set}
    \item{\code{clustsz}:}{Object of class \code{"list"} describing
      the number of instances for each cluster for each k tested}
    \item{\code{components}:}{Object of class \code{"vector"} storing
      the class membership of each protein for each k tested.}
    \item{\code{fcol}:}{Object of class \code{"character"} showing
      the feature column name in the corresponding \code{MSnSet}
      where the protein set information is stored.}
  }
} 

\section{Methods}{
  \describe{
    \item{plot}{Plots the kmeans clustering results.}
    \item{show}{Shows the object.}
  }
}


\author{
  Lisa M Breckels <lms79@cam.ac.uk>
}


\examples{
  showClass("ClustDist")
  
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
  
  ## plot clusters for first 'ClustDist' object 
  ## in the 'ClustDistList'
  plot(dd[[1]], xx)
  
  ## plot distances for all protein sets 
  plot(dd)
}

\keyword{classes}