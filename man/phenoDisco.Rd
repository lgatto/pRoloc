\name{phenoDisco}
\alias{phenoDisco}
\title{Runs the \code{phenoDisco} algorithm.}
\usage{
  phenoDisco(MSnSetObject, fcol = "markers", times = 10,
    GS = 6, allIter = FALSE, p = 0.05, r = 1, seed)
}
\arguments{
  \item{MSnSetObject}{An instance of class \code{MSnSet}.}

  \item{fcol}{A \code{character} indicating the organellar
  markers column name in feature meta-data. Default is
  \code{markers}.}

  \item{times}{Number of runs of tracking. Default is 10.}

  \item{GS}{Group size, i.e how many proteins make a group.
  Default is 6.}

  \item{allIter}{\code{logical}, defining if predictions
  for all iterations should be saved. Default is
  \code{FALSE}.}

  \item{p}{Significance level for outlier detection.
  Default is 0.05.}

  \item{r}{Correlation coefficent for Jaccard's Index.
  Default is 1.}

  \item{seed}{An optional \code{numeric} of length 1
  specifing the random number generator seed to be used.}
}
\value{
  An instance of class \code{MSnSet} containg the
  \code{phenoDisco} predictions.
}
\description{
  \code{phenoDisco} is a semi-supervised iterative approach
  to predict protein localisation and detect new protein
  clusters. It is based on the works of Yin et al. (2008).
}
\details{
  TODO - add description of the algorithm.
}
\author{
  Lisa M. Simpson <lms79@cam.ac.uk>
}
\references{
  Yin Z, Zhou X, Bakal C, Li F, Sun Y, Perrimon N, Wong ST.
  Using iterative cluster merging with improved gap
  statistics to perform online phenotype discovery in the
  context of high-throughput RNAi screens. BMC
  Bioinformatics. 2008 Jun 5;9:264. PubMed PMID: 18534020;
  PubMed Central PMCID: PMC2443381.
}

