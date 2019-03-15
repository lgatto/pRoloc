setClass("GenRegRes",
         representation(algorithm = "character",
                        hyperparameters = "list",
                        design = "numeric",
                        log = "list",
                        seed = "integer",
                        results = "matrix",
                        f1Matrices = "list", ## these are F1 matrices
                        cmMatrices = "list", ## these are continengy mat
                        testPartitions = "list",
                        datasize = "list"))

setMethod("show",
          signature(object = "GenRegRes"),
          function(object) {
            cat("Object of class \"",class(object),"\"\n",sep="")
            cat("Algorithm:", object@algorithm, "\n")
            cat("Hyper-parameters:\n")
            for (i in 1:length(object@hyperparameters)) {
              cat(" ", names(object@hyperparameters)[i],": ",
                  paste(object@hyperparameters[[i]], collapse = " "),
                  "\n", sep = "")
            }
            cat("Design:\n")
            cat(" Replication: ",
                object@design["times"], " x ",
                object@design["xval"], "-fold X-validation\n",
                sep = "")

            cat(" Partitioning: ",
                object@design["test.size"], "/",
                1-object@design["test.size"], " (test/train)\n",
                sep = "")

            cat("Results\n")
            res <- object@results
            cat(" macro F1:\n")
            print(summary(res[, "F1"]))
            for (i in 2:ncol(res)) {
              cat(" best ", colnames(res)[i],": ", sep = "")
              cat(paste(unique(res[, i])), collapse = " ", "\n")
            }
            if ("warnings" %in% names(object@log)) {
              cat("Use getWarnings() to see warnings.\n")
            }
            invisible(NULL)
          })

setMethod("getWarnings", "GenRegRes",
          function(object) {
            w <- object@log$warnings
            if (is.null(w)) {
              message("No warnings")
              invisible(w)
            } else {
              return(w)
            }
          })

setMethod("getSeed", "GenRegRes", function(object) object@seed)

setMethod("getF1Scores", "GenRegRes", function(object) object@results)

setMethod("f1Count", "GenRegRes",
          function(object, t) {
            f1tab <- getF1Scores(object)
            .f1 <- colnames(f1tab) == "F1"
            f1 <- f1tab[, .f1]
            if (missing(t))
              t <- max(f1, na.rm = TRUE)
            p <- f1tab[, !.f1, drop = FALSE]
            if (ncol(p) == 1) {
              res <- table(p[f1 >= t, ])
            } else {
              ## if ncol(p) != 1, then 2
              res <- tapply(f1,
                            list(factor(p[, 1]), factor(p[, 2])),
                            function(x) sum(x >= t))
            }
            return(res)
          })

setMethod("getParams", "GenRegRes",
          function(object) {
            res <- object@results
            best <- which.max(res[, "F1"])
            return(res[best, -1])
          })

setMethod("getOtherParams", "GenRegRes",
          function(object) {
            object@hyperparameters$other
          })

setMethod("plot", c("GenRegRes", "missing"),
          function(x, y, ...) {
            d <- data.frame(x@results)
            cn <- colnames(d)
            if (ncol(d) == 2) {
              colnames(d) <- c("F1", "p")
              p <- bwplot(F1 ~ factor(p),
                          data = d, xlab = cn[2])
            } else {
              colnames(d) <- c("F1", "p1", "p2")
              lp1 <- length(unique(d$p1))
              lp2 <- length(unique(d$p2))
              ifelse(lp1 < lp2,
                     p <- bwplot(F1 ~ factor(p2) | factor(p1), data = d, xlab = cn[3]),
                     p <- bwplot(F1 ~ factor(p1) | factor(p2), data = d, xlab = cn[2]))
            }
            p <- update(p,
                        main = x@algorithm,
                        panel = function(...) {
                          panel.grid(h = -1, v = 0)
                          panel.bwplot(...)
                        })
            p
          })


setMethod("levelPlot", "GenRegRes",
          function(object, fun = mean, ...) {
            x <- summariseMatList(object@f1Matrices, fun)
            labs <- names(dimnames(x))
            p <- levelplot(x, ylab = labs[2], xlab = labs[1],
                           main = object@algorithm,
                           ...)
            p
          })
