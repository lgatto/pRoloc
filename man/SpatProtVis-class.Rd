\name{SpatProtVis-class}
\Rdversion{1.1}
\docType{class}
\alias{SpatProtVis-class}
\alias{class:SpatProtVis}
\alias{SpatProtVis}
\alias{plot,SpatProtVis,missing-method}
\alias{show,SpatProtVis-method}

\title{Class \code{SpatProtVis}}

\description{

  A class for spatial proteomics visualisation, that upon instantiation,
  pre-computes all defined visualisations. Objects can be created with
  the \code{SpatProtVis} constructor and visualised with the \code{plot}
  method.

  The class is essentially a wrapper around several calls to
  \code{\link{plot2D}} that stores the dimensionality reduction
  outputs, and is likely to be updated in the future.
  
}

\usage{SpatProtVis(x, methods, dims, methargs, ...)}

\arguments{
  
  \item{x}{An instance of class \code{\linkS4class{MSnSet}} to visualise.}

  \item{methods}{Dimensionality reduction methods to be used to
    visualise the data. Must be contained in \code{\link{plot2Dmethods}}
    (except \code{"scree"}). See \code{\link{plot2D}} for details.}
  
  \item{dims}{A list of numerics defining dimensions used for
    plotting. Default are \code{1} and \code{2}. If provided, the length
    of this list must be identical to the length of \code{methods}. }
    
  \item{methargs}{A list of additional arguments to be passed for each
    visualisation method. If provided, the length of this list must be
    identical to the length of \code{methods}. }
  
  \item{\dots}{Additional arguments. Currently ignored.}
}

\section{Slots}{
  \describe{
   
    \item{\code{vismats}:}{A \code{"list"} of matrices containing the
      feature projections in 2 dimensions. }
    
    \item{\code{data}:}{The original spatial proteomics data stored as
      an \code{"MSnSet"}. }

    \item{\code{methargs}:}{A \code{"list"} of additional plotting
    arguments. }

    \item{\code{objname}:}{A \code{"character"} defining how to name the
      dataset. By default, this is set using the variable name used at
      object creation. }
    
  }
}

\section{Methods}{
  \describe{
    
    \item{\code{plot}:}{Generates the figures for the respective
      \code{methods} and additional arguments defined in the
      constructor. If used in an interactive session, the user is
      prompted to press 'Return' before new figures are displayed. }

    \item{\code{show}:}{A simple textual summary of the object.}
  
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\seealso{  
  The data for the individual visualisations is created by
  \code{\link{plot2D}}.
}

\examples{
library("pRolocdata")
data(dunkley2006)
## in the interest of time, don't use t-SNE
m <- c("PCA", "MDS", "kpca")
vis <- SpatProtVis(dunkley2006, methods = m)
vis
plot(vis)
plot(vis, legend = "topleft")

## Multiple PCA plots but different PCs
dims <- list(c(1, 2), c(3, 4))
vis <- SpatProtVis(dunkley2006, methods = c("PCA", "PCA"), dims = dims)
plot(vis)
}

\keyword{classes}
