
## References:
## [1] T. Burger, Y. Kessentini, T. Paquet. "Dempster-Shafer theory
## based rejection strategy for handwritten word recognition",
## The Eleventh International Conference on Document Analysis and Recognition
## (ICDAR 2011), Beijing, China, September 2011.

diffscores <- function(object, scols = "svm.scores") {
  ## See Diff in [1]
  ## Non continous values
  ## To be plotted and select threshold
  stopifnot(inherits(object, "MSnSet"))
  sidx <- grep(scols, fvarLabels(object))
  if (length(sidx) == 0)
    stop("No score 'svm.scores' columns found.")
  if (length(sidx) == 1) {
    msg <- c("\nOnly 1 ", fvarLabels(object)[sidx] , " column found.\n",
             "Are you sure you specified 'scores = \"all\"' when ",
             "running you classification?")
    stop(msg)
  }
  smat <- fData(object)[, sidx]
  ans <- apply(smat, 1, function(x) {
    .prob <- sort(x, decreasing = TRUE)[1:2]
    (.prob[1] - .prob[2])/.prob[1]    
  })
  if (missing(diffcol))
    diffcol <- sub("scores", "diff", scols)
  fData(x)[, diffcol] <- ans
  if (validObject(x))
    return(x)
}


viction <- function(object, scols = "svm.scores") {
  ## See Viction in [1]
  ## Phat is consonant mass function (masse de croyance)
  ## i is 1 .. Nunber of classes
  ## Phat_i = 1 * (P_i - P_(i+1)),
  ##   where P_i are ranked from largest to lowest prob 
  ## pl(A) = mass plaisibility = 1
  ## bel(A) = belief = Sum_i Phat_i
  ## where A: 1
  ##          1, 2
  ##          1, 2, ..., N
  ## Viction = Sum_A (1 - bel_A)   
  stopifnot(inherits(object, "MSnSet"))
  sidx <- grep(scols, fvarLabels(object))
  if (length(sidx) == 0)
    stop("No score 'svm.scores' columns found.")
  if (length(sidx) == 1) {
    msg <- c("\nOnly 1 ", fvarLabels(object)[sidx] , " column found.\n",
             "Are you sure you specified 'scores = \"all\"' when ",
             "running you classification?")
    stop(msg)
  }
  ## stub
}


difftorand <- function(object, scols = "svm.scores") {
  ## difference to random
  ## If we have N = length(sidx) classes, random assignment 
  ## can be approximiated by P_rand = 1/N [*] 
  ## (1) One could set as threshold P_i >= (n * P_rand).
  ## (2) Transform each score P_i2 = P_i - P_rand.
  ##     Positives score are more likely than random while
  ##     Negatives are less likely.
  ## (3) Calculate ratios P_i/P_rand
  ## May be all this is rather aimed as describing the
  ## scores. Especially a combination of (1) + (2) could
  ## be represented as barplots for each protein.
  ## 
  ## [*] This does not take into account class imbalance,
  ## which might screw these estimates up. May be use
  ## mean group_scores?   
  stopifnot(inherits(object, "MSnSet"))
  sidx <- grep(scols, fvarLabels(object))
  if (length(sidx) == 0)
    stop("No score 'svm.scores' columns found.")
  if (length(sidx) == 1) {
    msg <- c("\nOnly 1 ", fvarLabels(object)[sidx] , " column found.\n",
             "Are you sure you specified 'scores = \"all\"' when ",
             "running you classification?")
    stop(msg)
  }
  ## stub
}

