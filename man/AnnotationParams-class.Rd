\name{AnnotationParams-class}
\Rdversion{1.1}
\docType{class}
\alias{class:AnnotationParams}
\alias{AnnotationParams}
\alias{AnnotationParams-class}
\alias{show,AnnotationParams-method}
\alias{AnnotationParams}
\alias{setAnnotationParams}
\alias{getAnnotationParams}

\alias{dunkley2006params}
\alias{andy2011params}

\title{Class \code{"AnnotationParams"}}

\description{
  Class to store annotation parameters to automatically query a Biomart
  server, retrieve relevant annotation for a set of features of
  interest using, for example \code{\link{getGOFromFeatures}} and
  \code{\link{makeGoSet}}.
  
}


\section{Objects from the Class}{
  
  Objects can be created and set with the \code{setAnnotationParams}
  function. Object are created by calling without any arguments
  \code{setAnnotationParams()}, which will open an interactive
  interface. Depending on the value of \code{"many.graphics"} option, a
  graphical of a text-based menu will open (the text interface can be
  forced by setting the \code{graphics} argument to \code{FALSE}:
  \code{setAnnotationParams(graphics = FALSE)}). The menu will allow to
  select the species of interest first and the type of features (ENSEMBL
  gene identifier, Entrez id, ...) second.

  The species that are available are those for which ENSEMBL data is
  available in Biomart and have a set of attributes of interest
  available. The compatible identifiers for downstream queries are then
  automatically filtered and displayed for user selection.

  It is also possible to pass a parameter \code{inputs}, a character
  vector of length 2 containing a pattern uniquely matching the species
  of interest (in position 1) and a patterns uniquely matching the
  feature types (in position 2). If the matches are not unique, an error
  will be thrown.
  
  A new instance of the \code{AnnotationParams} will be created
  to enable easy and automatic query of the \code{Mart} instance. The
  instance is invisibly returned and stored in a global variable in the
  \pkg{pRoloc} package's private environment for automatic retrieval. If
  a variable containing an \code{AnnotationParams} instance is already
  available, it can be set globally by passing it as argument to the
  \code{setAnnotationParams} function. Globally set
  \code{AnnotationParams} instances can be accessed with the
  \code{getAnnotationParams} function.

  See the \code{pRoloc-theta} vignette for details.
  
}


\section{Slots}{
  \describe{

    \item{\code{mart}:}{Object of class \code{"Mart"} from the
      \pkg{biomaRt} package. }
    
    \item{\code{martname}:}{Object of class \code{"character"} with the
      name of the \code{mart} instance. }
    
    \item{\code{dataset}:}{Object of class \code{"character"} with the
      data set of the \code{mart} instance. }
    
    \item{\code{filter}:}{Object of class \code{"character"} with the
      filter to be used when querying the \code{mart} instance. }
    
    \item{\code{date}:}{Object of class \code{"character"} indicating
      when the current instance was created. }

   
    \item{\code{biomaRtVersion}:}{Object of class \code{"character"}
      with the \pkg{biomaRt} version used to create the
      \code{AnnotationParams} instance. }
   
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      with the version of the \code{AnnotationParams} class of the
      current instance. }
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "AnnotationParams")}: to
      display objects. }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
  
}

\seealso{
  
  \code{\link{getGOFromFeatures}}, \code{\link{makeGoSet}} and the
  \code{pRoloc-theta} vignette.

}

\examples{
data(andy2011params)
andy2011params
data(dunkley2006params)
dunkley2006params

try(setAnnotationParams(inputs = c("nomatch1", "nomatch2")))
setAnnotationParams(inputs = c("Homo sapiens",
                        "UniProtKB/Swiss-Prot ID"))
getAnnotationParams()
}

\keyword{classes}
