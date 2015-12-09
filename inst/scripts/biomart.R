## Last check Tue Dec  8 20:56:27 2015
## Also changed to use single mart instances

library("biomaRt")

## get MartInterface inferface
source("../../R/MartInterface.R")

mil <- MartInstanceList(list(ensembl = MartInstance(name = "ENSEMBL_MART_ENSEMBL",
                                                    host="www.ensembl.org",
                                                    path="/biomart/martservice"),
                             plants = MartInstance(name = "plants_mart",
                                                   host="plants.ensembl.org",
                                                   path="/biomart/martservice"),
                             fungi = MartInstance(name = "fungal_mart",
                                                  host="fungi.ensembl.org",
                                                  path="/biomart/martservice"),     
                             metazoa = MartInstance("metazoa_mart",
                                                    host="metazoa.ensembl.org",
                                                    path="/biomart/martservice")))
## filter out datasets that don't have all required attributes
mil <- lapply(mil, filterAttrs)

saveRDS(martTab, file="../extdata/mil.rds", compress = "xz")
