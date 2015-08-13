\name{SpatProtVis-class}
\Rdversion{1.1}
\docType{class}
\alias{SpatProtVis-class}
\alias{SpatProtVis}
\alias{plot,SpatProtVis,missing-method}
\alias{show,SpatProtVis-method}

\title{Class \code{SpatProtVis}}

\description{

A class for spatial proteomics visualisation, that upon instantiation,
pre-computes all available visualisations. Objects can be created with
the \code{SpatProtVis} constructor and visualised with the \code{plot}
method.

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

\section{Slots}{
  \describe{
    \item{\code{plot}:}{}
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
