getBestParams <- function(x) {
  if (all(is.na(x))) {
    msg <- paste("[pRoloc:::getBestParams] Only NA's in F1 matrix.\n",
                 " Try to use better suited parameters or check the marker class sizes with 'testMarkers(object)'.")
    warning(msg)
    x[1] <- 0 ## hack 
  } else if (any(is.na(x))) {
      warning("[pRoloc:::getBestParams] Found NA's in F1 matrix.")
  }
  k <- arrayInd(which( x == max(x, na.rm = TRUE) ),
                dim(x),
                useNames = TRUE)
  k <- apply(k, 1,
             function(i) as.numeric(c(colnames(x)[i["col"]],
                                      rownames(x)[i["row"]]))) 
  rownames(k) <- rev(names(dimnames(x)))  
  return(k)
}


makeF1matrix <- function(params) {
  ## The first item in the params list
  ## is used to create the columns,
  ## the second, if present, the rows.
  ## If there is only one item, the
  ## final matrix has one row.
  if (length(params) == 1) {
    .nrow <- 1
    .rnames <- "1"
  } else {
    .nrow <- length(params[[2]])
    .rnames <- params[[2]]
  }
  .ncol <- length(params[[1]])
  ans <- matrix(0, nrow = .nrow, ncol = .ncol)
  rownames(ans) <- .rnames
  colnames(ans) <- params[[1]]
  if (length(params) == 1) {
    names(dimnames(ans)) <- c("",
                              names(params)[1])
  } else {
    names(dimnames(ans)) <- c(names(params)[2],
                              names(params)[1])
  }
  return(ans)
}


checkNumbers <- function(x, tag, params) {
  sel <- is.nan(x)
  if (any(sel)) {      
    x[sel] <- 0
    if (!missing(tag)) {    
      new_warning <- paste0("NaNs found in '", tag, "' with hyperparameters ",
                            paste(names(params), params, sep = ":", collapse = " "),
                            ".")
      .warnings <- get(".warnings", envir = parent.frame())
      .warnings <- c(.warnings, new_warning)
      assign(".warnings", .warnings, envir = parent.frame())
    }
  }
  return(x)
}

subsetAsDataFrame <- function(object, fcol,
                              train = TRUE,
                              keepColNames = FALSE,
                              unknown = "unknown") {
    nms <- sampleNames(object)
    d <- data.frame(exprs(object), markers = fData(object)[, fcol])
    d.train <- d[d$markers != unknown,]
    d.train$markers <- factor(d.train$markers)
    d.test <- d[d$markers == unknown,]
    d.test$markers <- factor(d.test$markers)
    if (keepColNames)
        colnames(d) <- c(nms, fcol) 
    ifelse(train, 
           return(d.train),
           return(d.test))         
}

summariseMatList <- function(matList, fun = mean, ...) {
  stopifnot(length(unique(sapply(matList, ncol))) == 1)
  stopifnot(length(unique(sapply(matList, nrow))) == 1)
  ans <- apply(array(do.call(cbind, matList),
                     dim = c(nrow(matList[[1]]),
                       ncol(matList[[1]]),
                       length(matList))),
               c(1:2),
               fun, na.rm=TRUE, ...) ## assuming there is an na.rm argument in fun
  dimnames(ans) <- dimnames(matList[[1]])  
  ans
}

makePartList <- function(n, x) {
  .mkList2 <- function(i, x) {
    .res <- vector("list", length = x)
    names(.res) <- paste0("xval", 1:x)
    return(.res)
  }
  res <- lapply(seq_len(n), .mkList2, x)
  names(res) <- paste0("n", seq_len(n))
  return(res)
}
