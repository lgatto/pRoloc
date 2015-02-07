library(biomaRt)
library(pRoloc)

## Last check Fri Feb  6 21:42:58 2015

marts.of.interest <- c("ensembl","fungi_mart_25","metazoa_mart_25","plants_mart_25")

## > head(listMarts())
##               biomart                             version
## 1             ensembl        ENSEMBL GENES 78 (SANGER UK)
## 2                 snp    ENSEMBL VARIATION 78 (SANGER UK)
## 3 functional_genomics   ENSEMBL REGULATION 78 (SANGER UK)
## 4                vega                VEGA 58  (SANGER UK)
## 5       fungi_mart_25           ENSEMBL FUNGI 25 (EBI UK)
## 6 fungi_variations_25 ENSEMBL FUNGI VARIATION 25 (EBI UK)

martList <- lapply(marts.of.interest,
                   function(x) {
                     mart <- useMart(x)
                     ds <- listDatasets(mart)
                     cbind(ds, mart = x)
                   })

martTab <- do.call("rbind", martList)

## head(martTab)
##                   dataset                             description     version    mart
## 1  oanatinus_gene_ensembl  Ornithorhynchus anatinus genes (OANA5)       OANA5 ensembl
## 2   tguttata_gene_ensembl Taeniopygia guttata genes (taeGut3.2.4) taeGut3.2.4 ensembl
## 3 cporcellus_gene_ensembl         Cavia porcellus genes (cavPor3)     cavPor3 ensembl
## 4 gaculeatus_gene_ensembl  Gasterosteus aculeatus genes (BROADS1)     BROADS1 ensembl
## 5  lafricana_gene_ensembl      Loxodonta africana genes (loxAfr3)     loxAfr3 ensembl
## 6 mlucifugus_gene_ensembl        Myotis lucifugus genes (myoLuc2)     myoLuc2 ensembl


(n <- nrow(martTab)) ## 215
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
##                             name                       description
## 1                ensembl_gene_id                   Ensembl Gene ID
## 2          ensembl_transcript_id             Ensembl Transcript ID
## 3             ensembl_peptide_id                Ensembl Protein ID
## 4 canonical_transcript_stable_id Canonical transcript stable ID(s)
## 5                    description                       Description
## 6                chromosome_name                   Chromosome Name

head(filterList[[1]])
##                 name         description
## 1    chromosome_name     Chromosome name
## 2              start     Gene Start (bp)
## 3                end       Gene End (bp)
## 4             strand              Strand
## 5 chromosomal_region  Chromosome Regions
## 6      with_wikigene with WikiGene ID(s)

checkAttr0 <- function(x, attr0 = pRoloc2:::getAttributesOfInterest0()) 
    attr0 %in% x$name


checkAttrX <- function(x, attrX = pRoloc2:::getAttributesOfInterestX()) 
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
colnames(tmp) <- c(pRoloc2:::attributesOfInterest0,
                   sapply(pRoloc2:::attributesOfInterestX, "[", 1))

sel <- which(apply(tmp, 1, all))

attrList <- attrList[sel]
attr(attrList, "date") <- date()

filterList <- filterList[sel]
attr(filterList, "date") <- date()

martTab <- martTab[sel,]
attr(martTab, "date") <- date()

saveRDS(attrList, file="../extdata/attrList.rds")
saveRDS(filterList, file="../extdata/filterList.rds")
saveRDS(martTab, file="../extdata/martTab.rds")

