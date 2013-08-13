\name{phenoDisco}
\alias{phenoDisco}
\title{Runs the \code{phenoDisco} algorithm.}
\usage{
  phenoDisco(object, fcol = "markers", times = 100,
    GS = 10, allIter = FALSE, p = 0.05, seed,
    verbose = TRUE)
}
\arguments{
  \item{object}{An instance of class \code{MSnSet}.}

  \item{fcol}{A \code{character} indicating the organellar
  markers column name in feature meta-data. Default is
  \code{markers}.}

  \item{times}{Number of runs of tracking. Default is 100.}

  \item{GS}{Group size, i.e how many proteins make a group.
  Default is 10 (the minimum group size is 4).}

  \item{allIter}{\code{logical}, defining if predictions
  for all iterations should be saved. Default is
  \code{FALSE}.}

  \item{p}{Significance level for outlier detection.
  Default is 0.05.}

  \item{seed}{An optional \code{numeric} of length 1
  specifing the random number generator seed to be used.}

  \item{verbose}{Logical, indicating if messages are to be
  printed out during execution of the algorithm.}
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
  The algorithm performs a phenotype discovery analysis as
  described in Breckels et al. Using this approach one can
  identify putative subcellular groupings in organelle
  proteomics experiments for more comprehensive validation
  in an unbiased fashion. The method is based on the work
  of Yin et al. and used iterated rounds of Gaussian
  Mixture Modelling using the Expectation Maximisation
  algorithm combined with a non-parametric outlier
  detection test to identify new phenotype clusters.

  One requires 2 or more classes to be labelled in the data
  and at a very minimum of 6 markers per class to run the
  algorithm. The function will check and remove feature
  with missing values using the \code{\link{filterNA}}
  method.

  Important: Prior to version 1.1.2 the row order in the
  output was different from the row order in the input.
  This has now been fixed and row ordering is now the same
  in both input and output objects.
}
\examples{
\dontrun{
library(pRolocdata)
data(tan2009r1)
pdres <- phenoDisco(tan2009r1, fcol = "PLSDA")
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

  Breckels LM, Gatto L, Christoforou A, Groen AJ, Lilley KS
  and Trotter MWB. The Effect of Organelle Discovery upon
  Sub-Cellular Protein Localisation. J Proteomics. In
  Press.
}

