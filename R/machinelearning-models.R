
setMethod("consensusPrediction",
          "MSnSet",
          function(object,
                   method = c("filteredMajority", "majority", "strict"),
                   fcols = c("knn", "ksvm",
                     "naiveBayes",
                     "nnet", "plsda", "rf",
                     "svm"),
                   fun = mean) {
            ## fcols defines the predicters/scores
            ## to be used
            method <- match.arg(method)
            ans <- switch(method,
                          filteredMajority = majorityConsensus(object, fcols, filtered = TRUE, fun),
                          majority = majorityConsensus(object, fcols, filtered = FALSE),
                          strict = strictConsensus(object, fcols))
            return(ans)
          })


majorityConsensus <- function(object, fcols, filtered, fun) {
  preds <- sapply(fcols,
                  function(fcol) {
                    if (fcol %in% fvarLabels(object)) {
                      preds <- as.character(fData(object)[, fcol])
                      if (filtered) {
                        scrs <- fData(object)[, paste0(fcol, ".scores")]
                        k <- apply(scrs, 2, fun)
                        preds[scrs < k] <- "unknown"
                      }
                      return(preds)
                    }
                  })
  preds <- preds[!sapply(preds, is.null)]
  
  if (length(preds) == 0) {
    warning("No prediction results found using ",
            paste(fcols, collapse = ", "), ".")
  } else {
    stopifnot(length(unique(sapply(preds, length))) == 1)             
    preds <- do.call(cbind, preds)    
    preds <- apply(preds, 1, function(p) {
      tb <- table(p)
      x <- names(tb)[tb == max(tb)]
      ## If there are multiple majority predictions,
      ## 'unknown' are removed
      if (length(x) > 1)
        x <- x[x != "unknown"]
      nms <- paste(x, collapse = "|")
      nb <- length(x)
      if (filtered) {
        ans <- data.frame("filteredMajorityPreds" = nms,
                          "nbFilteredMajPreds" = nb)
        } else {
          ans <- data.frame("majorityPreds" = nms,
                            "nbMajPreds" = nb)
        }
      return(ans)
    })    
    preds <- do.call(rbind, preds)
    fData(object) <- cbind(fData(object), preds)
  }
  if (validObject(object))
    return(object)            
}


strictConsensus <- function(object, fcols) {
  preds <- sapply(fcols,
                  function(fcol) {
                    if (fcol %in% fvarLabels(object)) {
                      preds <- as.character(fData(object)[, fcol])
                      return(preds)
                    }
                  })
  preds <- preds[!sapply(preds, is.null)]
  if (length(preds) == 0) {
    warning("No prediction results found using ",
            paste(fcols, collapse = ", "), ".")
  } else {
    stopifnot(length(unique(sapply(preds, length))) == 1)             
    preds <- do.call(cbind, preds)
    preds <- apply(preds, 1,
                   function(p) {
                     up <- unique(p)
                     ifelse(length(up) == 1,
                            up,
                            "unknown")
                   })
    fData(object)$strictPreds <- preds
  }
  if (validObject(object))
    return(object)            
}



######################################
## tmp/not exported 

## a named list of GenRegRes models
getF1DataFrame <- function(x) {
  if (is.null(names(x)))
    names(x) <- sapply(x, function(.x) .x@algorithm)  
  xx <- c()
  for (n in names(x)) {
    tmp <- data.frame(F1 = x[[n]]@results[,"F1"],
                      model = n)
    xx <- rbind(xx, tmp)
  }
  return(xx)
}
