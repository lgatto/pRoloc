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
  computing \eqn{Chi^2}-values, the \eqn{\frac{sum (x_i–x_p)^2}{x_p}},
  where \eqn{x_i} is the normalised value of the peptide in fraction i
  and \eqn{x_p} is the value of the marker (Wiese et al 2007). The
  protein \eqn{Chi^2} is then computed as the median of the peptide
  \eqn{Chi^2} values. Peptides and proteins with similar profiles to the
  markers will have small \eqn{Chi^2} values. 

  The \code{chi2} methods implement this idea and compute such Chi^2
  values for sets of proteins.
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "matrix", y = "matrix")}}{
      Compute \code{nrow(x)} times \code{nrow(y)} \eqn{Chi^2} values,
      for each \code{x}, \code{y} feature pair.  
    }
    
    \item{\code{signature(x = "matrix", y = "numeric")}}{
      Computes \code{nrow(x)} \eqn{Chi^2} values, for all the \eqn{(x_i,
      y)} pairs.
    }
    
    \item{\code{signature(x = "numeric", y = "matrix")}}{
      Computes \code{nrow(y)} \eqn{Chi^2} values, for all the \eqn{(x,
      y_i)} pairs.

    }
    
    \item{\code{signature(x = "numeric", y = "numeric")}}{
      Computes the \eqn{Chi^2} value for the \eqn{(x, y)} pairs
    }
  }
}

\references{
  Andersen, J. S., Wilkinson, C. J., Mayor, T., Mortensen, P. et al.,
  Proteomic characterization of the human centrosome by protein
  correlation profiling. Nature 2003, 426, 570–574.

  Wiese, S., Gronemeyer, T., Ofman, R., Kunze, M. et al., Proteomics
  characterization of mouse kidney peroxisomes by tandem mass
  spectrometry and protein correlation profiling. Mol. Cell. Proteomics
  2007, 6, 2045–2057. 

}

\keyword{methods}

