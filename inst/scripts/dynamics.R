## msnl <- commonFeatureNames(list(n2014 = nightingale2014, n2015 = nightingale2015))
## d <- pRoloc:::pdist(msnl)

## d <- pRoloc:::pdist(yeastgg2, "condition")
## n <- 5
## o <- order(d, decreasing = TRUE)
## d[o[1:n]]

## yggfoi <- FeaturesOfInterest(description = "top dists", fnames = names(d[o])[1:n])

## p <- plot2D(yeastgg2[, yeastgg2$condition == "Glu"])
## highlightOnPlot(yeastgg2[, yeastgg2$condition == "Glu"],
##                 yggfoi, pch = 19, cex = 2)
## text(p[foi(yggfoi), ], labels = 1:n, col = "white")


## p <- plot2D(yeastgg2[, yeastgg2$condition == "Gly"])
## highlightOnPlot(yeastgg2[, yeastgg2$condition == "Gly"],
##                 yggfoi, pch = 19, cex = 2)
## text(p[foi(yggfoi), ], labels = 1:n, col = "white")


## yyg1 <- yeastgg2[, yeastgg2$replicate == 1]
## yyg2 <- yeastgg2[, yeastgg2$replicate == 2]
## d1 <- pRoloc:::pdist(yyg1, "condition")
## d2 <- pRoloc:::pdist(yyg2, "condition")

## plot((d1+d2)/2, log2(d1/d2))
## abline(h = 0)
