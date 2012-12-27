## machinelearning-framework
setGeneric("getWarnings", function(object) standardGeneric("getWarnings"))
setGeneric("getSeed", function(object) standardGeneric("getSeed"))
setGeneric("getF1Scores", function(object) standardGeneric("getF1Scores"))
setGeneric("levelPlot", function(object, ...) standardGeneric("levelPlot"))
setGeneric("getParams", function(object) standardGeneric("getParams"))
setGeneric("getRegularisedParams", function(object) getParams(object))
setGeneric("getRegularizedParams", function(object) getParams(object))

setGeneric("chi2",function(x,y,...) standardGeneric("chi2"))

setGeneric("consensusPrediction",function(object,...) standardGeneric("consensusPrediction"))
