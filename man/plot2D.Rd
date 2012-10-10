\name{plot2D}
\alias{plot2D}
\alias{plotDist}
\alias{addLegend}

\title{
  Plot organelle assignment data and results.
}

\description{
  Generate 2 dimansional or feature distribution plots to illustrate 
  localistation clusters. In \code{plot2D}, rows containing \code{NA}
  values are removed prior to dimention reduction.
}

\usage{
plot2D(object, fcol = "markers", fpch, unknown = "unknown", dims = 1:2,
alpha, score = 1, method = c("PCA", "MDS"), axsSwitch = FALSE, mirrorX =
FALSE, mirrorY = FALSE, col, pch, ...)

plotDist(object, markers, mcol = "steelblue", pcol = "grey90", alpha =
0.3, fractions, ...) 

addLegend(object, fcol = "markers", where = "other", ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An instance of class \code{MSnSet}. }
  \item{fcol}{ }
  \item{fpch}{ }
  \item{unknown}{ }
  \item{dims}{ }
  \item{alpha}{ }
  \item{score}{ }
  \item{method}{ }
  \item{axsSwitch}{ }
  \item{mirrorX}{ }
  \item{mirrorY}{ }
  \item{col}{ }
  \item{pch}{ }
  \item{where}{ (addLegend only) }
  \item{markers}{ }
  \item{mcol}{ }
  \item{pcol}{ }
  \item{fractions}{ }
  \item{\dots}{ }
}

\details{
  TODO  
}

\value{
  Used for their side effects of generating plots and adding legends.
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
library(pRolocdata)
data(dunkley2006)
plot2D(dunkley2006, fcol = NULL)
plot2D(dunkley2006, fcol = "markers")
addLegend(dunkley2006,
          fcol = "markers",
          where = "topright",
          cex = 0.5, bty = "n", ncol = 3)
title(main = "plot2D example")
data(tan2009r1)
j <- which(fData(tan2009r1)$markers == "mitochondrion")
i <- which(fData(tan2009r1)$PLSDA == "mitochondrion")
plotDist(tan2009r1[i, ],
         markers = featureNames(tan2009r1)[j],
         main = "Mitochondrion")
}