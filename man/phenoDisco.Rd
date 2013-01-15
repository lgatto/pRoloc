\name{phenoDisco}
\alias{phenoDisco}
\title{Runs the \code{phenoDisco} algorithm.}
\usage{
  phenoDisco(object, fcol = "markers", times = 50, GS = 10,
    allIter = FALSE, p = 0.05, r = 1, seed)
}
\arguments{
  \item{object}{An instance of class \code{MSnSet}.}

  \item{fcol}{A \code{character} indicating the organellar
  markers column name in feature meta-data. Default is
  \code{markers}.}

  \item{times}{Number of runs of tracking. Default is 50.}

  \item{GS}{Group size, i.e how many proteins make a group.
  Default is 10.}

  \item{allIter}{\code{logical}, defining if predictions
  for all iterations should be saved. Default is
  \code{FALSE}.}

  \item{p}{Significance level for outlier detection.
  Default is 0.05.}

  \item{r}{Correlation coefficent for Jaccard's Index.
  Default is 1 (currently no other value is supported).}

  \item{seed}{An optional \code{numeric} of length 1
  specifing the random number generator seed to be used.}
}
\value{
  An instance of class \code{MSnSet} containg the
  \code{phenoDisco} predictions.
}
\description{
  \code{phenoDisco} is a semi-supervised iterative approach
  to detect new protein clusters.
}
\details{
  The algorithm performs a phenotype discovery analysis
  which is based on the works of Yin et al. (2008). Using
  this appoach one can identify putative subcellular
  groupings in organelle proteomics experiments for more
  comprehensive validation in an unbiased fashion. The
  method uses iterated rounds of Gaussian Mixture Modelling
  using the Expectation Maximisation algorithm combined
  with a non-parametric outlier detection test to identify
  new phenotype clusters.
}
\examples{
\dontrun{
library(pRolocdata)
data(tan2009r1)
pdres <- phenoDisco(tan2009r1, markers = "PLSDA")
getPredictions(pdres, fcol = "pd", scol = NULL)
plot2D(pdres, fcol = "pd")
}
}
\author{
  Lisa M. Breckels <lms79@cam.ac.uk>
}
\references{
  Yin Z, Zhou X, Bakal C, Li F, Sun Y, Perrimon N, Wong ST.
  Using iterative cluster merging with improved gap
  statistics to perform online phenotype discovery in the
  context of high-throughput RNAi screens. BMC
  Bioinformatics. 2008 Jun 5;9:264. PubMed PMID: 18534020;
  PubMed Central PMCID: PMC2443381.
}

