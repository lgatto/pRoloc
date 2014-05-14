\name{chi2-methods}
\docType{methods}
\alias{chi2}
\alias{chi2-methods}
\alias{chi2,matrix,matrix-method}
\alias{chi2,matrix,numeric-method}
\alias{chi2,numeric,matrix-method}
\alias{chi2,numeric,numeric-method}

\title{The PCP 'chi square' method}

\description{
  In the original protein correlation profiling (PCP), Andersen et
  al. use the peptide normalised profiles along gradient fractions and
  compared them with the reference profiles (or set of profiles) by
  computing \eqn{Chi^2} values,
  \eqn{\frac{\sum (x_i - x_p)^2}{x_p}}{sum(x_i - x_p)^2 / x_p},
  where \eqn{x_i} is the normalised value of the peptide in fraction i
  and \eqn{x_p} is the value of the marker (from Wiese et al., 2007). The
  protein \eqn{Chi^2} is then computed as the median of the peptide
  \eqn{Chi^2} values. Peptides and proteins with similar profiles to the
  markers will have small \eqn{Chi^2} values. 

  The \code{chi2} methods implement this idea and compute such Chi^2
  values for sets of proteins. 
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "matrix", y = "matrix", method =
	"character", fun = "NULL", na.rm = "logical")}}{
      Compute \code{nrow(x)} times \code{nrow(y)} \eqn{Chi^2} values,
      for each \code{x}, \code{y} feature pair. Method is one of
      \code{"Andersen2003"} or \code{"Wiese2007"}; the former (default)
      computed the \eqn{Chi^2} as \code{sum(y-x)^2/length(x)}, while the
      latter uses \code{sum((y-x)^2/x)}. \code{na.rm} defines if missing
      values (\code{NA} and \code{NaN}) should be removed prior to
      summation. \code{fun} defines how to summarise the \eqn{Chi^2}
      values; default, \code{NULL}, does not combine the \eqn{Chi^2}
      values. 
    }
    
    \item{\code{signature(x = "matrix", y = "numeric", method =
	"character", na.rm = "logical")}}{
      Computes \code{nrow(x)} \eqn{Chi^2} values, for all the \eqn{(x_i,
	y)} pairs. See above for the other arguments.
    }
    
    \item{\code{signature(x = "numeric", y = "matrix", method =
	"character", na.rm = "logical")}}{
      Computes \code{nrow(y)} \eqn{Chi^2} values, for all the \eqn{(x,
	y_i)} pairs.  See above for the other arguments.
    }
    
    \item{\code{signature(x = "numeric", y = "numeric", method =
	"character", na.rm = "logical")}}{
      Computes the \eqn{Chi^2} value for the \eqn{(x, y)} pairs.  See
      above for the other arguments. 
    }
  }
}

\references{
  Andersen, J. S., Wilkinson, C. J., Mayor, T., Mortensen, P. et al., 
  Proteomic characterization of the human centrosome by protein 
  correlation profiling. Nature 2003, 426, 570 - 574.

  Wiese, S., Gronemeyer, T., Ofman, R., Kunze, M. et al., Proteomics
  characterization of mouse kidney peroxisomes by tandem mass
  spectrometry and protein correlation profiling. Mol. Cell. Proteomics
  2007, 6, 2045 - 2057. 

}

\seealso{ \code{\link{empPvalues}} }

\examples{
mrk <- rnorm(6)
prot <- matrix(rnorm(60), ncol = 6)
chi2(mrk, prot, method = "Andersen2003")
chi2(mrk, prot, method = "Wiese2007")

pepmark <- matrix(rnorm(18), ncol = 6)
pepprot <- matrix(rnorm(60), ncol = 6)
chi2(pepmark, pepprot)
chi2(pepmark, pepprot, fun = sum)
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\keyword{methods}

