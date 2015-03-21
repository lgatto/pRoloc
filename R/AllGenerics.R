## machinelearning-framework
setGeneric("getWarnings", function(object) standardGeneric("getWarnings"))
setGeneric("getSeed", function(object) standardGeneric("getSeed"))
setGeneric("getF1Scores", function(object) standardGeneric("getF1Scores"))
setGeneric("levelPlot", function(object, ...) standardGeneric("levelPlot"))
setGeneric("getParams", function(object, ...) standardGeneric("getParams"))

setGeneric("getRegularisedParams",
           function(object) {
               .Deprecated("getParams", "pRoloc")
               getParams(object)
           })

setGeneric("getRegularizedParams",
           function(object) {
               .Deprecated("getParams", "pRoloc")
               getParams(object)
           })

setGeneric("getOtherParams", function(object, ...) standardGeneric("getOtherParams"))

setGeneric("f1Count", function(object, ...) standardGeneric("f1Count"))

setGeneric("chi2", function(x,y,...) standardGeneric("chi2"))

## clustering
setGeneric("kmeansClustering", function(object, params, ...) standardGeneric("kmeansClustering"))
setGeneric("kmeansOptimisation", function(object, ...) standardGeneric("kmeansOptimisation"))
setGeneric("kmeansOptimization", function(object, ...) kmeansOptimisation(object, ...))

## distances
setGeneric("nndist", function(object, query, ...) standardGeneric("nndist"))

