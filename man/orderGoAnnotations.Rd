% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/goannotations.R
\name{orderGoAnnotations}
\alias{orderGoAnnotations}
\title{Orders annotation information}
\usage{
orderGoAnnotations(
  object,
  fcol = "GOAnnotations",
  k = 1:5,
  n = 5,
  p = 1/3,
  verbose = TRUE,
  seed
)
}
\arguments{
\item{object}{An instance of class \code{MSnSet}.}

\item{fcol}{The name of the annotations matrix. Default is
\code{GOAnnotations}.}

\item{k}{The number of clusters to test. Default is \code{k = 1:5}}

\item{n}{The minimum number of proteins per component cluster.}

\item{p}{The normalisation factor, per \code{k} tested}

\item{verbose}{A \code{logical} indicating if a progress bar should
be displayed. Default is \code{TRUE}.}

\item{seed}{An optional random number generation seed.}
}
\value{
An updated \code{MSnSet} containing the newly ordered 
\code{fcol} matrix.
}
\description{
For a given matrix of annotation information, 
this function returns the information ordered according to 
the best fit with the data.
}
\details{
As there are typically many protein/annotation sets that may fit the data
we order protein sets by best fit i.e. cluster tightness, by
computing the mean normalised Euclidean distance for all instances 
per protein set. 

For each protein set i.e. proteins that have been labelled
with a specified term/information criteria, we find the best 
\code{k} cluster components for the set (the default is to 
test\code{k = 1:5}) according to the minimum mean normalised 
pairwise Euclidean distance over all component clusters. 
(Note: when testing \code{k} if any components are found to 
have less than \code{n} proteins these components are not
included and \code{k} is reduced by 1). 

Each component cluster is normalised by \code{N^p} (where 
\code{N} is the total number of proteins per component, 
and \code{p} is the power). Hueristally, \code{p = 1/3} 
and normalising by \code{N^1/3} has been found the optimum 
normalisation factor. 

Candidates in the matrix are ordered according to lowest 
mean normalised pairwise Euclidean distance as we expect 
high density, tight clusters to have the smallest mean 
normalised distance.

This function is a wrapper for running \code{clustDist},
\code{getNormDist}, see the "Annotating spatial proteomics data"
vignette for more details.
}
\seealso{
\code{addGoAnnotations} and example therein.
}
\author{
Lisa M Breckels
}
