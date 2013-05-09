setClass("ClustRegRes",
         representation(algorithm = "character"))

setMethod("show",
          signature(object = "ClustRegRes"),
          function(object) {
            cat("Object of class \"",class(object),"\"\n",sep="")
            cat("Algorithm:", object@algorithm, "\n")
            invisible(NULL)            
          })
