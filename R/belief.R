difftorand <- function(object, scols = "svm.scores") {
  ## difference (or ratio) to random
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
    stop("No '", scols , "' scores columns found.")
  if (length(sidx) == 1) {
    msg <- c("\nOnly 1 ", fvarLabels(object)[sidx] , " column found.\n",
             "Are you sure you specified 'scores = \"all\"' when ",
             "running you classification?")
    stop(msg)
  }
  smat <- fData(object)[, sidx]
  prand <- 1/ncol(smat)
  as.matrix(smat / prand)
}

.diffscore <- function(p) {
  if (length(p) < 2)
    stop("must have a least 2 probabilities to compute the diff scores")
  p <- sort(p, decreasing = TRUE)[1:2]
  (p[1] - p[2])/p[1]    
}

diffscores <- function(object, scols = "svm.scores") {
  ## See Diff in [1]
  ## Non continous values
  ## To be plotted and select threshold
  stopifnot(inherits(object, "MSnSet"))
  sidx <- grep(scols, fvarLabels(object))
  if (length(sidx) == 0)
    stop("No '", scols , "' scores columns found.")
  if (length(sidx) == 1) {
    msg <- c("\nOnly 1 ", fvarLabels(object)[sidx] , " column found.\n",
             "Are you sure you specified 'scores = \"all\"' when ",
             "running you classification?")
    stop(msg)
  }
  smat <- fData(object)[, sidx]
  ans <- apply(smat, 1, .diffscore)
  ## if (missing(diffcol))
  ##   diffcol <- sub("scores", "diff", scols)
  ## fData(x)[, diffcol] <- ans
  ## if (validObject(x))
  ##   return(x)
  return(ans)
}


phat <- function(p) {
  ## Calculates the consonant mass
  ## function a sorted vector of
  ## probabilities
  stopifnot(all.equal(sum(p), 1))
  stopifnot(all(p >= 0 & p <= 1))
  p <- sort(p)
  n <- length(p)
  ans <- numeric(n) 
  p <- c(sort(p, decreasing = TRUE), 0)  
  for (i in 1:n)
    ans[i] <- i * (p[i] - p[i+1])
  ans
}

powerset <- function(n) {
  if (length(n) > 1)
    n <- length(n)
  ans <- sapply(1:n,
                function(i)
                combn(n, i, simplify = FALSE))
  unlist(ans, recursive = FALSE)
}


pl <- function(Phat, A, B) {
  ## pl(A) = \sum_{(B \cap A) \not= \emptyset} m(B), \forall A \subseteq \Omega
  if (missing(A))
    A <- powerset(Phat)
  if (missing(B))
    B <- sapply(1:length(Phat), function(i) c(1:i))
  idx <- sapply(A, function(a)
                sapply(B, function(b)
                       ifelse(length(intersect(b, a) > 0),
                              TRUE, FALSE)))  
  apply(idx, 2, function(i) sum(Phat[i]))
}

bel <- function(Phat = NULL, A, B) {
  ## pl(A) = \sum_{(B \subseteq A), B \not= \emptyset} m(B), \forall A \subseteq \Omega
  if (missing(A))
    A <- powerset(Phat)
  if (missing(B))
    B <- sapply(1:length(Phat), function(i) c(1:i))
  idx <- sapply(A, function(a)
                sapply(B, function(b)
                       ifelse(all(is.element(b, a)),
                              TRUE, FALSE)))
  apply(idx, 2, function(i) sum(Phat[i]))
}


.viction <- function(P) {
  Phat <- phat(P)
  A <- powerset(Phat)
  B <- sapply(1:length(Phat), function(i) c(1:i))
  .pl <- pl(Phat, A, B)
  .bel <- bel(Phat, A, B)
  sum(.pl - .bel)
}


viction <- function(object, scols = "svm.scores") {
  ## See Viction in [1]
  ## Phat is consonant mass function
  ## i is 1 .. Number of classes
  ## P = sort(scores)
  ## Phat_i = 1 * (P_i - P_(i+1)),
  ## pl(A) = mass plaisibility = 1
  ## bel(A) = belief = Sum_i Phat_i
  ## where A: 1
  ##          1, 2
  ##          1, 2, ..., N
  ## Viction = Sum_A (1 - bel_A)   
  stopifnot(inherits(object, "MSnSet"))
  sidx <- grep(scols, fvarLabels(object))
  if (length(sidx) == 0)
    stop("No '", scols , "' scores columns found.")
  if (length(sidx) == 1) {
    msg <- c("\nOnly 1 ", fvarLabels(object)[sidx] , " column found.\n",
             "Are you sure you specified 'scores = \"all\"' when ",
             "running you classification?")
    stop(msg)
  }
  smat <- fData(object)[, sidx]
  apply(smat, 1, .viction)
}


## References:
## [1] T. Burger, Y. Kessentini, T. Paquet. "Dempster-Shafer theory
## based rejection strategy for handwritten word recognition",
## The Eleventh International Conference on Document Analysis and Recognition
## (ICDAR 2011), Beijing, China, September 2011.
