\name{nndist-methods}
\docType{methods}
\alias{nndist-methods}
\alias{nndist,matrix,matrix-method}
\alias{nndist,matrix,missing-method}
\alias{nndist,MSnSet,missing-method}
\alias{nndist}

\title{ Nearest neighbour distances }

\description{
  Methods computing the nearest neighbour indices and distances for
  \code{matrix} and \code{MSnSet} instances.   
}

\section{Methods}{

  \describe{

    \item{\code{signature(object = "matrix", k = "numeric", dist =
	"character", ...)}}{ Calculates indices and distances to the
	\code{k} (default is 3) nearest neighbours of each feature (row)
	in the input matrix \code{object}. The distance \code{dist} can
	be either of \code{"euclidean"} or
	\code{"mahalanobis"}. Additional parameters can be passed to the
	internal function \code{FNN::get.knn}. Output is a matrix with
	\code{2 * k} columns and \code{nrow(object)} rows. }

    \item{\code{signature(object = "MSnSet", k = "numeric", dist =
	"character", ...)}}{ As above, but for an \code{MSnSet}
	input. The indices and distances to the \code{k} nearest
	neighbours are added to the object's feature metadata.  }
    
    \item{\code{signature(object = "matrix", query = "matrix", k =
      "numeric", ...)}}{ If two \code{matrix} instances are provided as
      input, the \code{k} (default is 3) indices and distances of the
      nearest neighbours of \code{query} in \code{object} are returned
      as a matrix of dimensions \code{2 * k} by
      \code{nrow(query)}. Additional parameters are passed to
      \code{FNN::get.knnx}. Only euclidean distance is available. }

  }
}

\keyword{methods}

\examples{
library("pRolocdata")
data(dunkley2006)

## Using a matrix as input
m <- exprs(dunkley2006)
m[1:4, 1:3]
head(nndist(m, k = 5))
tail(nndist(m[1:100, ], k = 2, dist = "mahalanobis"))

## Same as above for MSnSet
d <- nndist(dunkley2006, k = 5)
head(fData(d))

d <- nndist(dunkley2006[1:100, ], k = 2, dist = "mahalanobis")
tail(fData(d))

## Using a query
nndist(m[1:100, ], m[101:110, ], k = 2)
}
