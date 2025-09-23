getFiltersOfInterest <- function()
   .Defunct()
    # c("ensembl_gene_id",
    #   "ensembl_transcript_id",
    #   "ensembl_peptide_id",
    #   "ensembl_exon_id",
    #   "wikigene_id",
    #   "wikigene_name",
    #   "embl",
    #   "entrezgene",
    #   "hgnc_id",
    #   "hgnc_symbol",
    #   "protein_id",
    #   "refseq_dna",
    #   "refseq_mrna",
    #   "refseq_peptide",
    #   "refseq_genomic",
    #   "ottt",
    #   "ottg",
    #   "ucsc",
    #   "uniprotsptrembl",
    #   "uniprotswissprot",
    #   "uniprot_gn",
    #   "unigene",
    #   "geneindex",
    #   "tair_locus",
    #   "tair_locus_model",
    #   "tair_symbol",
    #   "flybase_annotation_id",
    #   "flybase_gene_id",
    #   "flybase_transcript_id",
    #   "flybase_translation_id",
    #   "flybasename_gene",
    #   "flybasename_transcript",
    #   "flybasename_translation",
    #   "sgd",
    #   "sgd_gene",
    #   "sgd_transcript")

getGONamespace <- function()
    .Defunct()
    # c("biological_process", "cellular_component",
    #   "molecular_function")

getAttributesOfInterest0 <- function()
    .Defunct()
    # c("ensembl_gene_id", "ensembl_peptide_id",
    #   "ensembl_transcript_id",
    #   "cdna", "peptide", "coding",
    #   "gene_exon", "gene_exon_intron",
    #   "transcript_exon_intron",
    #   "go_linkage_type")

getAttributesOfInterestX <- function()
    .Defunct()
    # list(c("go_id", "go_accession", "goslim_goa_accession"),
    #      c("namespace_1003", "go_namespace_1003"),
    #      c("name_1006", "go_name_1006"))

chooseAttributesX <-
  function(p) {
    .Defunct()
    # .attrX <- getAttributesOfInterestX()
    # sapply(.attrX, function(.attr) {
    #   sel <- .attr %in% biomaRt::listAttributes(p@mart)[, 1]
    #   .attr[sel][1] ## return first in case multiple matches
    # })
  }
