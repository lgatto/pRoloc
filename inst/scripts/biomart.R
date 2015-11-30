library(biomaRt)
library(pRoloc)

## Updated on "Mon Nov 30 22:06:33 2015"

## > biomaRt::listMarts()
##                biomart               version
## 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 82
## 2     ENSEMBL_MART_SNP  Ensembl Variation 82
## 3 ENSEMBL_MART_FUNCGEN Ensembl Regulation 82
## 4    ENSEMBL_MART_VEGA               Vega 62
## 5                pride        PRIDE (EBI UK)


## 1         ensembl ENSEMBL GENES 79 (SANGER UK)
## 5   fungi_mart_26    ENSEMBL FUNGI 26 (EBI UK)
## 7 metazoa_mart_26  ENSEMBL METAZOA 26 (EBI UK)
## 9  plants_mart_26   ENSEMBL PLANTS 26 (EBI UK)

marts.of.interest <- "ENSEMBL_MART_ENSEMBL"

martList <- lapply(marts.of.interest,
                   function(x) {
                     mart <- useMart(x)
                     ds <- listDatasets(mart)
                     cbind(ds, mart = x)
                   })

martTab <- do.call("rbind", martList)

## head(martTab)
##                          dataset                                description
## 1         oanatinus_gene_ensembl     Ornithorhynchus anatinus genes (OANA5)
## 2        cporcellus_gene_ensembl            Cavia porcellus genes (cavPor3)
## 3        gaculeatus_gene_ensembl     Gasterosteus aculeatus genes (BROADS1)
## 4         lafricana_gene_ensembl         Loxodonta africana genes (loxAfr3)
## 5 itridecemlineatus_gene_ensembl Ictidomys tridecemlineatus genes (spetri2)
## 6        choffmanni_gene_ensembl        Choloepus hoffmanni genes (choHof1)
##   version                 mart
## 1   OANA5 ENSEMBL_MART_ENSEMBL
## 2 cavPor3 ENSEMBL_MART_ENSEMBL
## 3 BROADS1 ENSEMBL_MART_ENSEMBL
## 4 loxAfr3 ENSEMBL_MART_ENSEMBL
## 5 spetri2 ENSEMBL_MART_ENSEMBL
## 6 choHof1 ENSEMBL_MART_ENSEMBL

(n <- nrow(martTab)) ## 69
attrList <- filterList <- vector("list", length = n)
names(attrList) <- names(filterList) <- martTab$dataset[1:n]

pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
  setTxtProgressBar(pb, i)
  martname <- as.character(martTab[i,"mart"])
  dataset <- as.character(martTab[i,"dataset"])
  mart <- useMart(martname,dataset)
  attrList[[i]] <- listAttributes(mart)
  filterList[[i]] <- listFilters(mart)
}
close(pb)


head(attrList[[1]])
##                    name           description         page
## 1       ensembl_gene_id       Ensembl Gene ID feature_page
## 2 ensembl_transcript_id Ensembl Transcript ID feature_page
## 3    ensembl_peptide_id    Ensembl Protein ID feature_page
## 4       ensembl_exon_id       Ensembl Exon ID feature_page
## 5           description           Description feature_page
## 6       chromosome_name       Chromosome Name feature_page


head(filterList[[1]])
##                 name                                               description
## 1    chromosome_name                                           Chromosome name
## 2              start                                           Gene Start (bp)
## 3                end                                             Gene End (bp)
## 4             strand                                                    Strand
## 5 chromosomal_region Chromosome Regions (e.g 1:100:10000:-1,1:100000:200000:1)
## 6          with_hgnc                                           with HGNC ID(s)

checkAttr0 <- function(x, attr0 = pRoloc:::getAttributesOfInterest0()) 
    attr0 %in% x$name


checkAttrX <- function(x, attrX = pRoloc:::getAttributesOfInterestX()) 
    sapply(attrX,
           function(.attrX) {
               stopifnot(length(.attrX) == 2)
               .attrX[1] %in% x$name | .attrX[2] %in% x$name
           })


## selecting only those mart data sets that
## have all the attributes of interest
tmp <- t(sapply(attrList, checkAttr0))
tmp2 <- t(sapply(attrList, checkAttrX))
tmp <- cbind(tmp, tmp2)

## colnames(tmp) <- c(attributes.of.interest0,                   
##                    sapply(attributes.of.interestX, "[", 1))
colnames(tmp) <- c(pRoloc:::getAttributesOfInterest0(),
                   sapply(pRoloc:::getAttributesOfInterestX(), "[", 1))

sel <- which(apply(tmp, 1, all))

attrList <- attrList[sel]
attr(attrList, "date") <- date()

filterList <- filterList[sel]
attr(filterList, "date") <- date()

martTab <- martTab[sel,]
attr(martTab, "date") <- date()

saveRDS(attrList, file="../extdata/attrList.rds", compress = "xz")
saveRDS(filterList, file="../extdata/filterList.rds", compress = "xz")
saveRDS(martTab, file="../extdata/martTab.rds", compress = "xz")
