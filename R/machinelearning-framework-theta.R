setClass("ThetaRegRes",
         contains = "GenRegRes",
         slots = list(
             predictions = "list", 
             otherWeights = "list"))

setMethod("show",
          signature(object = "ThetaRegRes"),
          function(object) {
              cat("Object of class \"",class(object),"\"\n",sep="")
              cat("Algorithm:", object@algorithm, "\n")
              cat("Theta hyper-parameters:\n")
              hyper <- object@hyperparameters$theta
              w <- unique(as.vector(hyper))
              w <- w[order(w)]
              wx <- round(w, 2)
              cat(" weights:", w, "\n") 
              cat(" k:", object@hyperparameters$k, "\n")
              cat(" nrow:", nrow(object@hyperparameters$theta), "\n")
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
              resF1 <- object@results[, "F1"]     
              resTh <- object@results[, c(2:ncol(object@results))]
              cat(" macro F1:\n")
              print(summary(resF1))
              w <- as.factor(w)
              if (is.vector(resTh)) {
                  m <- data.frame(t(data.frame(resTh)))
              } else {
                  m <- data.frame(resTh)
              }            
              ll <- levels(w)
              foo <- function(z) factor(z, levels = ll)
              m <- colwise(foo)(m)
              foo <- function(z) table(z)
              t <- colwise(foo)(m)
              rownames(t) <- paste("weight:", wx, sep = "")
              cat(" best theta:\n")
              print(t)
              if ("warnings" %in% names(object@log)) {
                  cat("Use getWarnings() to see warnings.\n")
              }
              invisible(NULL)            
          })

## object  =   Object of class "ThetaRegRes"
## t        =  Numeric. If specified thetas with an f1 score > t are considered.
## For thetaOpt results, we look at results and first select to examine results 
## only with F1 > t, then of those, we look to see if any thetas are duplicated,
## if yes we count
setMethod("f1Count", "ThetaRegRes",
          function(object, t) {
            f1tab <- getF1Scores(object)
            .f1 <- f1tab[, "F1"]
            .weights <- f1tab[, 2:ncol(f1tab), drop = FALSE]
            if (missing(t))
              t <- max(.f1, na.rm = TRUE)
            .tokeep <- which(.f1 >= t)
            if (length(.tokeep)==1) {
              res <- t(data.frame(.weights[.tokeep, ]))
              res <- c(res, 1)
              res <- matrix(res, nrow = 1)
              colnames(res) <- c(colnames(.weights), "count")
            } else {
              df <- data.frame(.weights[.tokeep, ])
              res <- ddply(df, colnames(df), nrow)
              res <- as.matrix(res)
              colnames(res) <- c(colnames(.weights), "count")
              .order <- order(res[, "count"], decreasing = TRUE)
              res <- res[.order, ]
            }            
            return(res)
          })


## object   =  Object of class "ThetaRegRes"
## t        =  Numeric. If specified thetas with an f1 score > t are considered.   
setMethod("getParams", "ThetaRegRes",
          function(object, 
                   method = c("median", "mean", "max", "count"),
                   favourP = FALSE
                   ) {
            if (missing(method)) method <- "median"

            if (method == "max") {
              scores <- object@results[, 1]
              ind <- which(scores == max(scores))
              if (length(ind) > 1) {
                res <- object@results[ind, -1, drop = FALSE]
                if (favourP) {
                  message("More than one best theta, picking weights that favour primary. See plot() to examine best thetas")
                  rs <- rowSums(res)
                  ind <- which(rs == max(rs))
                  if (length(ind) > 1) {
                    ind <- sample(ind, size = 1)
                  }
                } else {
                  message("More than one best theta, picking weights at random from best. See plot() to examine best thetas")
                  ind <- sample(1:nrow(res), size = 1)
                }
                best <- res[ind, ]
              } else {
                best <- object@results[ind, -1]
              }
            }
            
            if (method == "median") {
              res <- object@results[, -1, drop = FALSE]
              best <- apply(res, 2, median)
            }
            
            if (method == "mean") {
              res <- object@results[, -1, drop = FALSE]
              best <- apply(res, 2, mean)
            }
            
            if (method == "count") {
              res <- f1Count(object, t = 0)
              idx <- ncol(res)
              cc <- res[, idx]
              ind <- which(cc == max(cc))
              res <- res[ind, -idx]
              if (class(res) != "numeric") {
                if (favourP) {
                  message("More than one best theta, picking weights that favour primary. See plot() to examine best thetas")
                  rs <- rowSums(res)
                  ind <- which(rs == max(rs))
                  if (length(ind) > 1) {
                    ind <- sample(ind, size = 1)
                  }
                } else {
                  message("More than one best theta, picking weights at random from best. See plot() to examine best thetas")
                  ind <- sample(1:nrow(res), size = 1)
                }
                best <- res[ind, ]
              } else {
                best <- res
              }
            }
          
            return(best)
          })


##      x =  Object of call "ThetaRegRes"
## labels =  Logical. Indicates whether to print the frequency 
##           of weight assigments per class
##      n =  A numeric that specifies the maximum size of the plotting symbol.
##           Default is 20.
##      t =  Numeric. If specified thetas with an f1 score > t are considered.   
setMethod("plot", c("ThetaRegRes", "missing"),
          function(x, y,
                   n = 20,
                   labels = FALSE,
                   axsSwitch = FALSE,
                   cols = getStockcol(),
                   wlab, cex = 14,
                   t) {
              if (nrow(x@hyperparameters$theta) == 1)
                  stop("Only one theta weight vector tested, no weights to plot")
              if (!is.logical(labels))
                  stop("labels must be of class logical")
              if (!is.numeric(n))
                  stop("n must be of class numeric")
              if (missing(t))
                  t <- min(x@results[, "F1"])
              f1 <- x@results[, "F1"]
              tokeep <- which(f1 >= t)
              best <- x@results[tokeep, c(2:ncol(x@results))]
              th <- unique(as.vector(x@hyperparameters$theta))
              w <- as.vector(best)
              if (is.vector(best)) { 
                  orgs <- names(best)
                  orgs <- factor(orgs, levels = rev(orgs))
                  df <- data.frame(orgs, best, row.names=c(1:length(best)))
                  colnames(df) <- c("classes", "w")
              } else {
                  orgs <- colnames(best)
                  classes <- factor(rep(orgs, each = nrow(best)), levels = rev(orgs))
                  df <- data.frame(classes, w)
              }
              df[,"w"] <- factor(w, levels = th)
              data <- ddply(df, c("classes", "w"), "nrow", .drop = FALSE)
            
              nth <- length(th)
              ncl <- ncol(best)
              if (axsSwitch) {
                  aes <- aes(x = classes, y = w, size = nrow)
              } else {
                  aes <- aes(x = w, y = classes, size = nrow)
              }  
              ssd <- scale_size_continuous(range = c(0, n)) ## Scaling of points
              points <- geom_point(aes(color = classes))
              theme.par <- theme(legend.position = "none", 
                                 axis.text = element_text(size = cex,
                                     colour = rgb(0, 0, 0)), ## Axis ticks and labels font size
                                 text = element_text(size = cex + 1), ## Axis labels font size
                                 axis.title.y = element_text(vjust = 0.8, colour = rgb(0, 0, 0)), 
                                 axis.title.x = element_text(vjust = -.5, colour = rgb(0, 0, 0))) 
              ## Adjust space between axis and titles with vjust
              ## set.ticks <- scale_y_continuous(breaks=th, limits = c(-0.2, 1.2))
              ## q <- ggplot(data, aes) + points + ssd + theme.par + set.ticks + xlab("Class") + ylab("Classifier weight")
              if (missing(wlab)) 
                  wlab <- sort(round(th, 3))
              if (axsSwitch) {
                  q <- ggplot(data, aes) + points + ssd + theme.par + xlab("Class") + ylab("Classifier weight")
                  q <- q + scale_y_discrete(labels = wlab)
              } else {
                  q <- ggplot(data, aes) + points + ssd + theme.par +
                      xlab("Classifier weight") + ylab("Class")
                  q <- q + scale_x_discrete(labels = wlab)
              }            
              if (labels) {
                  q <- q + geom_text(data = subset(data, nrow != 0), 
                                     aes(label = nrow), size = 4)
              } 
              ## q <- q + scale_fill_manual(values = cols[1:ncl])
              q <- q + scale_color_manual(values = rev(cols[1:ncl]))
              q 
          })        


plotThetas <- function(object, 
                       smooth = TRUE,
                       medians = TRUE,
                       proportion = TRUE,
                       colramp,
                       ...) {
    stopifnot(inherits(object, "ThetaRegRes"))
    if (nrow(object@hyperparameters$theta) == 1)
        stop("Only one theta weight vector tested, no thetas to plot")          
    th <- rowSums(object@hyperparameters$theta)
    ## matrix where each row is a possible theta combination and columns are runs
    data <- matrix(unlist(object@f1Matrices), nrow = length(th))
    ## use medians: take the median f1 score over all runs for each combination
    if (medians) {
        f1 <- apply(data, 1, function(x) median(x))
        r <- range(f1)
    } else {
        f1 <- as.vector(data)
        th <- rep(th, ncol(data))
        r <- range(f1)
    }
    if (proportion) {
        th <- th/max(th)
        ylab <- "Proportion of primary to auxiliary data"
        w <- which.max(f1)
        loc <- th[w] + .02
    } else {
        ylab <- "Row sum of thetas"
    }
    if (missing(colramp)) colramp <- c("white", blues9)
    if (smooth) smoothScatter(f1, th, xlab = "F1 scores", 
                              xlim = r, ylab = ylab, 
                              colramp = colorRampPalette(colramp), ...)
    else plot(f1, th, xlab = "F1 scores", xlim = r,
              ylab = ylab, ...)
    max <- max(f1)
    w <- which.max(f1)
    if (!proportion) {
        myth <- object@hyperparameters$theta
        minT <- min(rowSums(myth))
        maxT <- max(rowSums(myth))
        abline(h=minT, lty=3, lwd=1.5)
        abline(h=maxT, lty=3, lwd=1.5)
        text(x=r[1]+.03, y=round(minT+.2, 5), cex=.6, 
             labels=paste("min(sum of thetas) = ", minT))
        text(x=r[1]+.03, y=round(maxT-.2, 3), cex=.6, 
             labels=paste("max(sum of thetas) = ", maxT))
        loc <- th[w]+.18
    }
    abline(h=th[w], lty=3, lwd=1.5)
    text(y=loc, x=r[1]+.03, cex=.6, 
         labels=paste("max (F1 score) = ", round(max, 3)))
}


combineThetaRegRes <- function(object) {
    ##  object: a list of "ThetaRegRes" objects
    if (!class(object)=="list")
        stop("Object must be of class list")

    if (!all(sapply(object, inherits, "ThetaRegRed")))
        stop("All elements in list must be of class ThetaRegRes")
    .thetaSlots <- slotNames("ThetaRegRes")
    .slotChk <- sapply(object, function(z) slotNames(z) == .thetaSlots)    
    
    if (any(!.slotChk))
        stop("Slot names do not match ThetaRegRes slots")
    
    .xval <- sapply(object, function(z) z@design["xval"])
    .test.size <- sapply(object, function(z) z@design["test.size"])
    if (length(unique(.xval)) > 1 | length(unique(.test.size)) > 1)
        stop("Every object must be run with same xval and test.size")
    
    .ks <- sapply(object, function(z) z@hyperparameters$k)
    .kChk<- apply(.ks, 1, function(z) length(unique(z)))
    if (any(.kChk > 1)) stop("The same k must be used for every run to combine")
    
    .thetas <- lapply(object, function(z) z@hyperparameters$theta)
    .tf <- sapply(.thetas, function(z) identical(.thetas[[1]], z))
    if (any(!.tf)) stop ("The same thetas must be used for every run to combine")
    
    .fcol <- sapply(object, function(z) z@datasize$fcol)
    if (length(unique(.fcol)) > 1)
        stop("The same fcol must be used in every run to combine")
    
    
    .results <- do.call(rbind, lapply(object, function(z) z@results))
    .hyperparams <- object[[1]]@hyperparameters
    .timesID <- sapply(object, function(z) z@design["times"])
    .times <- sum(.timesID)
    .design <- c(.xval[1], .test.size[1], .times)
    names(.design) <- c("xval", "test.size", "times")
    .algorithm <- object[[1]]@algorithm
    .predictions <- do.call(c, lapply(object, function(z) z@predictions))
    .cmMatrices <- do.call(c, lapply(object, function(z) z@cmMatrices))
    .f1Matrices <- do.call(c, lapply(object, function(z) z@f1Matrices))
    .names <- paste("times", 1:.times, sep="")
    names(.f1Matrices) <- .names
    .otherWeights <- do.call(c, lapply(object, function(z) z@otherWeights))
    .ds <- object[[1]]@datasize
    .test.partitions <- vector("list", 2)
    names(.test.partitions) <- c("validation", "train")
    .validation <- lapply(object, function(z) z@testPartitions$validation)
    .validation <- do.call(c, .validation)
    names(.validation) <- paste("validation", 1:.times, sep ="")
    .train <- lapply(object, function(z) z@testPartitions$train)
    .train <- do.call(c, .train)
    names(.train) <- paste("train", 1:.times, sep ="")
    .test.partitions[[1]] <- .validation
    .test.partitions[[2]] <- .train  
    .log <- vector("list", 2)
    names(.log) <- c("warnings", "splits")
    .tf <- which(sapply(.otherWeights, is.null) == FALSE)
    if (length(.tf) > 0) {
        .log[[1]] <- paste0("For ", length(.otherWeights), 
                            " of the times there are multiple best thetas")}
    .splits <- lapply(object, function(z) z@log$splits)
    .splits <- sapply(1:length(.timesID), function(z) 
        replicate(.timesID[z], .splits[[z]], simplify = FALSE))
    .splits <- do.call(c, .splits)
    .log$splits <- .splits
    
    ans <- new("ThetaRegRes",
               algorithm = .algorithm,
               hyperparameters = .hyperparams,
               design = .design,
               log = .log,
               results = .results,
               f1Matrices = .f1Matrices,
               cmMatrices = .cmMatrices,
               testPartitions = .test.partitions, 
               datasize = .ds,
               predictions = .predictions,
               otherWeights = .otherWeights)
    return(ans)
}
