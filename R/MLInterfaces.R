setMethod("MLearn",
          c("formula", "MSnSet", "learnerSchema","numeric" ),          
          function( formula, data, .method, trainInd, ...) {
            data <- as(data,"ExpressionSet")
            thecall <- match.call()
            ans <- MLearn(formula, data, .method, trainInd, ...)
            ans@call <- thecall
            ans@learnerSchema <- .method
            return(ans)
          })


setMethod("MLearn",
          c("formula", "MSnSet", "learnerSchema", "xvalSpec" ),
          function( formula, data, .method, trainInd, ...) {
            thecall <- match.call()
            data <- as(data,"ExpressionSet")
            data <- MLInterfaces:::es2df(data, keep=as.character(as.list(formula)[[2]]))
            ans <- MLearn(formula, data, .method, trainInd, ...)
            ans@call <- thecall
            ans@learnerSchema <- .method
            return(ans)
          })
