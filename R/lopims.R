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
                  message("[LOPIMS 2] Synergising file ", i, " out of ", n)
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

refNormMedOfPepRatio <- function(x) {
  ## Polpitiya, A., Qian, W., Jaitly, N., Petyuk, V. et al., 
  ## DAnTE: a statistical tool for quantitative analysis of -omics data. 
  ## Bioinformatics 2008, 24, 1556–1558.
  ##
  ## See also Matzke et al. Proteomics DOI 10.1002/pmic.201200269
  ## Ratio all peptides to the peptide with the least amount of missing values;
  ## Protein abundance is the median of the scaled peptide abundances
  ref <- which.min(apply(x, 2, function(.x) sum(is.na(.x))))
  x <- as.matrix(x)
  .res <- apply(x, 1, function(.x) {
    .xref <- .x[ref]
    .x/.xref
  })
  res <- t(.res)
  apply(res, 2, median, na.rm = TRUE)
}  

refNormMeanOfNonNAPepSum <- function(x) {
  ## For each protein:
  ## (1) Find ref fraction with least amount of missing values
  ## (2) For each (frac_i, frac_ref), consider only peptides
  ##     without missing values, break ties using max(sum(int))
  ## (3) Calculate mean(sum(frac_i), sum(frac_ref))
  ref <- apply(x, 2, function(.x) sum(is.na(.x)))
  idx <- which(ref == min(ref)) ## indices if min na columns 
  if (length(idx) > 1) {
    sms <- colSums(x, na.rm = TRUE)
    ## index of max sum from the min na columns
    jdx <- which.max(sms[idx]) 
    ref <- idx[jdx]
  } else {
    ref <- which.min(ref)
  }
  x <- as.matrix(x)
  .res <- rep(-1, ncol(x))
  for (i in 1:ncol(x)) {
    sel_i <- !is.na(x[, i])
    sel_ref <- !is.na(x[, ref])
    sel <- sel_i & sel_ref
    .res[i] <- sum(x[sel, i])/sum(x[sel, ref])
  }
  .res
}  

.lopims4 <- function(xx) {
  ## remove rows with only NAs
  xx <- xx[!apply(exprs(xx), 1, function(x) all(is.na(x))), ]
  message("[LOPIMS 4] Calculate ratios to reference peptide")
  res <- by(exprs(xx),
            fData(xx)$protein.Accession.LOPIMS1,
            refNormMeanOfNonNAPepSum)  
  res2 <- Reduce(rbind, res)
  rownames(res2) <- names(res)
  xx2 <- combineFeatures(xx, fData(xx)$protein.Accession.LOPIMS1, fun = median)
  ## xx2 <- MSnbase:::nologging(xx2) 
  exprs(xx2) <- res2
  xx2@processingData@processing <-
    c(xx2@processingData@processing,
      paste0("Normalised intensities to reference fraction: ",
             date()))
  sampleNames(xx2) <- sampleNames(xx)
  if (validObject(xx2))
    return(xx2)
}



.lopims5 <- function(x, markerfile) {
  message("[LOPIMS 5] Adding markers")
  addMarkers(x, markerfile, verbose = TRUE)
}

##' The function processes MSe data using the \code{\link{synergise}}
##' function of the \code{\link{synapter}} package and combines resulting
##' \code{\link{Synapter}} instances into one \code{"\linkS4class{MSnSet}"}
##' and organelle marker data is added as a feature-level annotation variable.
##' 
##' The \code{LOPIMS} pipeline is composed of 5 steps:
##' \enumerate{
##' \item The HDMSe final peptide files are used to compute 
##' false discovery rates uppon all possible combinations of
##' HDMSe final peptides files and the best combination smaller
##' or equal to \code{mfdr} is chosen.
##' See \code{\link{estimateMasterFdr}} for details.
##' The corresponding master run is then created as descibed in
##' \code{\link{makeMaster}}.
##' 
##' \item Each MSe/pep3D pair is processed using the HDMSe master
##' file using \code{\link{synergise}}.
##'
##' \item The respective peptide-level \code{synergise} output objects
##' are converted and combined into an single \code{"\linkS4class{MSnSet}"}
##' instance.
##'
##' \item Protein-level quantitation is inferred as follows.
##' For each protein, a reference sample/fraction is chosen based on the
##' number of missing values (\code{NA}). If several samples have a same
##' minimal number of \code{NA}s, ties are broken using the sum of counts.
##' The peptide that do not display any missing values for each
##' (frac_{i}, frac_{ref}) pair are summed and the ratio is reported
##' (see pRoloc:::refNormMeanOfNonNAPepSum for details).
##'
##' \item The markers defined in the \code{markerfile} are collated as
##' feature meta-data in the \code{markers} variable.
##' See \code{\link{addMarkers}} for details.
##' }
##'
##' Intermediate \code{synergise} reports as well as resulting objects
##' are stored in a \code{LOPIMS_pipeline} directory.
##' For details, please refer to the \code{synapter} vignette and
##' reference papers.
##'
##' @references
##' Improving qualitative and quantitative performance for MSE-based label free proteomics.
##' N.J. Bond, P.V. Shliaha, K.S. Lilley and L. Gatto, Journal of Proteome Research, 2013 (in press).
##'
##' The Effects of Travelling Wave Ion Mobility Separation on Data Independent Acquisition in Proteomics Studies,
##' P.V. Shliaha, N.J. Bond, L. Gatto and K.S. Lilley, Journal of Proteome Research, 2013 (in press).
##' 
##' MSnbase-an R/Bioconductor package for isobaric tagged mass spectrometry data visualization,
##' processing and quantitation. L. Gatto and KS. Lilley. Bioinformatics. 2012 Jan 15;28(2):288-9.
##' doi: 10.1093/bioinformatics/btr645. Epub 2011 Nov 22. PubMed PMID: 22113085.
##' 
##' @title A complete LOPIMS pipeline
##' @param hdmsedir A \code{character} identifying the directory
##' containing the HDMSe final peptide files. Default is \code{HDMSe}.
##' @param msedir A \code{character} identifying the directory
##' containing the MSe final peptide files. Default is \code{MSe}.
##' @param pep3ddir A \code{character} identifying the directory
##' containing the MSe pep 3D files. Default is \code{pep3D}.
##' @param fastafile A \code{character} identifying the protein
##' fasta database. Default is to use the fasta file in the current
##' directory. If several such files exist, the function reports an
##' error.
##' @param markerfile A \code{character} identifying the marker
##' file (see details for format). Default is to use a \code{csv} file
##' starting with \code{marker} in the current directory.
##' If several such files exist, the function reports an error.
##' @param mfdr The master FDR value. Default is 0.025.
##' @param ... Additional paramters passed to \code{\link{synergise}}.
##' @return An instance of class \code{"\linkS4class{MSnSet}"} with protein
##' level quantitation and respective organelle markers.
##' @author Laurent Gatto
lopims <- function(hdmsedir = "HDMSE",
                   msedir = "MSE",
                   pep3ddir = "pep3D",
                   fastafile,
                   markerfile,
                   mfdr = 0.025,
                   ...) {
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
    markerfile <- dir(pattern = "marker.*\\.csv$")
    if (length(markerfile) != 1)
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
  lopims <- .lopims5(msn, markerfile)
  write.exprs(lopims,
              file = file.path(mainoutputdir, "lopims.txt"),
              fDataCols = fvarLabels(lopims))
  save(lopims, file = file.path(mainoutputdir, "lopims.rda"))
  SOFTWARE <- file.path(mainoutputdir, "SOFTWARE")
  cat(msg0, file = SOFTWARE)
  cat(msg1, file = SOFTWARE, append = TRUE)
  invisible(lopims)
}
