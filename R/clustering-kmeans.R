setGeneric("kmeansClustering",
           function(object, params, ...)
           standardGeneric("kmeansClustering"))

setMethod("kmeansClustering",
          signature(object = "MSnSet",
                    params = "missing"),
          function(object, ...) {
            cl <- kmeans(exprs(object), ...)
            fData(object)$kmeans <- cl$cluster
            ## possibly update processingData
            if (validObject(object))
              return(object)
          })

setMethod("kmeansClustering",
          signature(object = "matrix",
                    params = "missing"),
          function(object, ...) {
            kmeans(object, ...)
          })

setMethod("kmeansClustering",
          signature(object = "MSnSet",
                    params = "ClustRegRes"),
          function(object, params){
            ## extract params
            ## object <- kmeans(object, params)
            ## return(object)          
          })

setMethod("kmeansOptimisation",
          signature(object = "MSnSet"),
          function(object, ...) {
            ## optimise and create/return a ClustGenRes
          })
