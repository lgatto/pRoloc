getFiltersOfInterest <- function()
    c("ensembl_gene_id",
      "ensembl_transcript_id",
      "ensembl_peptide_id",
      "ensembl_exon_id",
      "wikigene_id",
      "wikigene_name",
      "embl",
      "entrezgene",
      "hgnc_id",
      "hgnc_symbol",
      "ipi",
      "protein_id",
      "refseq_dna",
      "refseq_mrna",
      "refseq_peptide",
      "refseq_genomic",
      "ottt",
      "ottg",
      "ucsc",
      "uniprot_sptrembl",
      "uniprot_swissprot",
      "uniprot_swissprot_accession",
      "unigene",
      "geneindex",
      "tair_locus",
      "tair_locus_model",
      "tair_symbol",
      "uniprot_genename",
      "flybase_annotation_id",
      "flybase_gene_id",
      "flybase_transcript_id",
      "flybase_translation_id",
      "flybasename_gene",
      "flybasename_transcript",
      "flybasename_translation",
      "sgd",
      "sgd_gene",
      "sgd_transcript")

getGONamespace <- function()
    c("biological_process", "cellular_component",
      "molecular_function")

## These functions return the attributes that are used as a criterion
## for ENSEMBL mart sets to be used in pRoloc. These are the
## attributes that are susceptible to be queried by pRoloc functions
## as \code{\link{getGOFromFeatures}} or \code{\link{getSeq}}.
##
## @title Get Biomart attributes of interest. Returns a
## \code{character} or \code{list} of \code{characters} (when multiple
## alternatives exists) with the attributes that the biomart data sets
## used in pRoloc all have.  @author Laurent Gatto @rdname
## getAttributesOfInterest

getAttributesOfInterest0 <- function()
    c("ensembl_gene_id", "ensembl_peptide_id",
      "ensembl_transcript_id",
      "cdna", "peptide", "coding",
      "gene_exon", "gene_exon_intron",
      "transcript_exon_intron",
      "go_linkage_type")

getAttributesOfInterestX <- function()
    list(c("go_id", "go_accession"),
         c("namespace_1003", "go_namespace_1003"),
         c("name_1006", "go_name_1006"))

chooseAttributesX <-
  function(p) {  
    .attrX <- getAttributesOfInterestX()
    sapply(.attrX, function(.attr) {
      sel <- .attr %in% biomaRt::listAttributes(p@mart)[, 1]
      .attr[sel][1] ## return first in case multiple matches
    })
  }
