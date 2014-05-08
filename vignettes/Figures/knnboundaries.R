library("pRoloc")
library("pRolocdata")
data(tan2009r1)

## Working in 2D

x <- dunkley2006[, 1:2]
exprs(x) <- plot2D(dunkley2006, plot = FALSE)
sampleNames(x) <- paste0("PC", 1:2)


## ----------------------

x <- knnClassification(x, k = 1, fcol = "pd.markers")
fData(x)$knn1 <- fData(x)$knn
x <- knnClassification(x, k = 10, fcol = "pd.markers")
fData(x)$knn10 <- fData(x)$knn
x <- knnClassification(x, k = 20, fcol = "pd.markers")
fData(x)$knn20 <- fData(x)$knn

## kp <- knnOptimisation(dunkley2006, fcol = "pd.markers")
## load("dunkley2006knnpar.rda")
## best k = 5

x <- knnClassification(x, k = 5, fcol = "pd.markers")
fData(x)$knn5 <- fData(x)$knn
fData(x)$knn <- NULL


## x <- exprs(dunkley2006)
## pca <- prcomp(x, center = TRUE, scale = TRUE)
## px <- pca$x ## is x %*% pca$rotation

## unprcomp <- function(px, pca) {
##     x <- px %*% t(pca$rotation)
##     x <- sweep(x, 2L, pca$scale, "*")
##     x <- sweep(x, 2L, pca$center, "+")
##     x
## }

## stopifnot(all.equal(unprcomp(px, pca), x))
## stopifnot(all.equal(px, predict(pca, newdata = x)))


pc1 <- seq(-6.5, 7, 0.1)
pc2 <- seq(-6, 5, 0.1)
gm <- as.matrix(expand.grid(pc1, pc2))
rownames(gm) <- paste0("G", 1:nrow(gm))
colnames(gm) <- sampleNames(x)

cl <- factor(fData(x)$pd.markers)
tr <- fData(x)$pd.markers != "unknown"


cl1 <- class::knn(exprs(x)[tr, ], exprs(x)[!tr, ], cl[tr], k = 1)
cl5<- class::knn(exprs(x)[tr, ], exprs(x)[!tr, ], cl[tr], k = 5)
cl10 <- class::knn(exprs(x)[tr, ], exprs(x)[!tr, ], cl[tr], k = 10)

gcl1 <- class::knn(exprs(x)[tr, ], gm, cl[tr], k = 1)
gcl5 <- class::knn(exprs(x)[tr, ], gm, cl[tr], k = 5)
gcl8 <- class::knn(exprs(x)[tr, ], gm, cl[tr], k = 8)
gcl10 <- class::knn(exprs(x)[tr, ], gm, cl[tr], k = 10)

cls <- paste0(getStockcol(), 60)[-(6:7)]
cls0 <- getStockcol()[-(6:7)]
setUnknowncol("#00000010")

.plot <- function(.cl, .gcl, ...) {
    plot(exprs(x), ...)
    points(exprs(x)[tr, 1], exprs(x)[tr, 2], pch = 19, col = cls0[cl[tr]])    
    points(exprs(x)[!tr, 1], exprs(x)[!tr, 2], pch = 19, col = cls[.cl])
    points(gm[, 1], gm[, 2], col = cls[.gcl], pch = "+", cex = 0.8)
}


## pdf(width = 7, height = 7)
## plot(exprs(x))
## points(exprs(x)[tr, 1], exprs(x)[tr, 2], pch = 19, col = cls[cl[tr]])    
## points(gm[, 1], gm[, 2], col = cls0[gcl1], pch = "+", cex = 0.8)
## ## points(gm[, 1], gm[, 2], col = cls0[gcl10], pch = 1, cex = 0.8)
## dev.off()

x1 <- c(0, 3); y1 <- c(3, 0.15)
x2 <- c(1.7, 4); y2 <- c(4.4, 1)

pdf("knnboundaries.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
.plot(NA, gcl1, main = "Markers and k = 1 classification boundaries.")
rect(x1, y1, x2, y2, lwd = 2)
## .plot(NA, gcl10, main = "Markers and k = 10 classification boundaries.")
.plot(NA, gcl8, main = "Markers and k = 8 classification boundaries.")
rect(x1, y1, x2, y2, lwd = 2)
## .plot(cl1, gcl5, main = "knn 1")
## .plot(cl10, gcl5, main = "knn 10")
dev.off()

plot2D(dunkley2006, fcol = "pd.markers", main = "markers", pch = 1)
points(pgm[, 1], pgm[, 2], col = cls[gcl10], pch = "+")

plot2D(dunkley2006, fcol = "knn1", main = "knn 1", pch = 1)
points(pgm[, 1], pgm[, 2], col = cls[gcl5], pch = "+")
plot2D(dunkley2006, fcol = "knn20", main = "knn 20", pch = 1)
points(pgm[, 1], pgm[, 2], col = cls[gcl5], pch = "+")
plot2D(dunkley2006, fcol = "knn5", main = "knn 5", pch = 1)
points(pgm[, 1], pgm[, 2], col = cls[gcl5], pch = "+")

