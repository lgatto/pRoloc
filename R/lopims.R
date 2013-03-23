.lopims1 <- function(hdmsedir, fastafile, mfdr = 0.025) {
  HDMSeFinalPeptideFiles <- dir(hdmsedir, full.names = TRUE)
  message("[LOPIMS 1] Estimate master FDR")
  cmb <- estimateMasterFdr(HDMSeFinalPeptideFiles,
                           fastafile,
                           masterFdr = mfdr,
                           verbose = TRUE)
  print(cmb)
  message("[LOPIMS 1] Make master run")
  master <- makeMaster(HDMSeFinalPeptideFiles[bestComb(cmb)],
                       verbose = TRUE)
  print(master)
  return(master)
}

.lopims2 <- function(msedir, pep3ddir, fasta,
                     master, outdir,
                     ...) {
  ## master file
  masterfile <- file.path(outdir, "master.rds")
  saveRDS(master, file = masterfile)
  ## mse files
  msefiles <- dir(msedir, full.names = TRUE)
  pep3dfiles <- dir(pep3ddir, full.names = TRUE)
  stopifnot(length(msefiles) == length(pep3dfiles))
  n <- length(msefiles)
  ## iterate of each file pairs
  res <- lapply(seq_along(msefiles),
                function(i) {
                  message("[LOPIMS 2] Synergising files ", i, " out of ", n)
                  .outputdir <- file.path(outdir, paste0("synapter_results_", i))
                  dir.create(.outputdir)
                  .input <- list(identpeptide = masterfile,
                                 quantpeptide = msefiles[i],
                                 quantpep3d = pep3dfiles[i],
                                 fasta = fasta)
                  synergise(filenames = .input,
                            outputdir = .outputdir,
                            master = TRUE,
                            ...)
                })
  return(res)
}

.lopims3 <- function(synlist) {
  message("[LOPIMS 3] Combine synapter results")
  ll <- lapply(synlist, as, "MSnSet")
  nms <- paste0("LOPIMS", seq_along(ll))
  ll <- lapply(seq_along(ll), function(i) {
    sampleNames(ll[[i]]) <- nms[i]
    ll[[i]] <- updateFvarLabels(ll[[i]], nms[i])
    ll[[i]]
  })      
  Reduce(combine, ll)
}

.lopims4 <- function(xx) {
  ## remove rows with only NAs
  xx <- xx[!apply(exprs(xx), 1, function(x) all(is.na(x))), ]
  ref <- which.min(apply(exprs(xx), 2, function(x) sum(is.na(x))))
  message("[LOPIMS 4] Calculate ratios to reference column ", sampleNames(xx)[ref])
  res <- by(exprs(xx), fData(xx)$protein.Accession.LOPIMS1,
            function(x) {
              x <- as.matrix(x)
              .res <- apply(x, 1, function(.x) {
                .xref <- .x[ref]
                .x <- .x[-ref]
                .x/.xref
              })
              res <- t(.res)
              apply(res, 2, mean, na.rm = TRUE)
            })
  
  res2 <- Reduce(rbind, res)
  rownames(res2) <- names(res)

  xx2 <- combineFeatures(xx, fData(xx)$protein.Accession.LOPIMS1, fun = mean)
  xx2 <- xx2[, -ref]
  xx2 <- MSnbase:::nologging(xx2)
  exprs(xx2) <- res2
  if (validObject(xx2))
    return(xx2)
}

.lopims5 <- function(x, markerfile) {
  message("[LOPIMS 5] Adding markers")
  addMarkers(x, markerfile, verbose = TRUE)
}

lopims <- function(hdmsedir = "HDMSE",
                   msedir = "MSE",
                   pep3ddir = "pep3D",
                   fastafile,
                   markerfile,
                   mfdr = 0.025,
                   ...) {
  suppressPackageStartupMessages(library("MSnbase"))
  suppressPackageStartupMessages(library("synapter"))
  suppressPackageStartupMessages(library("pRoloc"))

  msg0 <-
    paste("\n ------------------------------------------------------------\n",
          "This free open-source software implements academic research. \n",
          "If you use it, please cite:\n\n",          
          "L. Gatto and KS. Lilley. MSnbase - an R/Bioconductor\n",
          "package for isobaric tagged mass spectrometry data visualization,\n",
          "processing and quantitation. Bioinformatics 28, 288-289 (2012).\n\n",          
          "NJ. Bond, PV. Shliaha, KS. Lilley and L. Gatto\n", 
          "Improving qualitative and quantitative performance for MSE-based\n",
          "label free proteomics Journal of Proteome Research, in press (2013).\n\n",          
          "L. Gatto, LM. Breckels with contributions from T. Burger and\n", 
          "S. Wieczorek. pRoloc: A unifying bioinformatics framework for \n",
          "spatial proteomics. R package (2013).\n",
          "------------------------------------------------------------\n\n")
  message(msg0)
  
  version <- "0.1.0"
  released <- "2013/03/22"  
  msg1 <- paste0("LOPIMS pipline version ", version, " (", released, "), using\n",
                 "synapter ", packageVersion("synapter"),
                 ", MSnbase ", packageVersion("MSnbase"), 
                 " and pRoloc ", packageVersion("pRoloc"), ".\n\n") 
  message(msg1)
  
  if (missing(fastafile)) {
    fastafile <- dir(pattern = "\\.fasta$")
    if (length(fastafile) != 1)
      stop("Require exactly one fasta file in directory.")
  }
  if (missing(markerfile)) {
    markerfile <- dir(pattern = "marker.*\\.fasta$")
    if (length(fastafile) != 1)
      stop("Require exactly one markers csv file in directory.")
  }
  stopifnot(file.exists(hdmsedir))
  stopifnot(file.exists(msedir))
  stopifnot(file.exists(pep3ddir))
  mainoutputdir <- paste0("LOPIMS_pipeline_", gsub(" ", "_", date()))
  dir.create(mainoutputdir)
  message("[LOPIMS  ] Results will be saved to ", mainoutputdir)  
  
  master <- .lopims1(hdmsedir, fastafile, mfdr)
  res <- .lopims2(msedir, pep3ddir, fastafile, master,
                  mainoutputdir, ...)
  msn <- .lopims3(res)
  msn <- .lopims4(msn)
  msn <- .lopims5(msn, markerfile)
  write.exprs(msn,
              file = file.path(mainoutputdir, "msn.txt"),
              fDataCols = fvarLabels(msn))
  save(msn, file = file.path(mainoutputdir, "msn.rda"))
  SOFTWARE <- file.path(mainoutputdir, "SOFTWARE")
  cat(msg0, file = SOFTWARE)
  cat(msg1, file = SOFTWARE, append = TRUE)
  invisible(msn)
}
