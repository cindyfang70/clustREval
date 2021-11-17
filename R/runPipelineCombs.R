library(pipeComp)
library(scuttle)
library(labelled)

#' runPipelineCombs
#' 
#' A function that runs user-defined clustering pipelines by taking all possible combinations of user-provided pipeline hyperparameters.
#' User-provided hyperparameters must match available options. 
#' 
#' 
#' @param sce SingleCellExperiment object with raw counts in counts slot.
#' @param outputPrefix String indicating the prefix of the filepath that results should save to
#' @param doubletmethod Character vector of doublet detection methods. Default: none
#' @param filt Character vector of filtering methods. Default: filt.default
#' @param norm Character vector of normalization methods Default: norm.seurat
#' @param sel Character vector of selection methods. Default: sel.vst
#' @param selnb Numeric vector of number of genes to select. Default: 2000
#' @param dr Character vector of dimensionality reduction methods. Default: seurat.pca
#' @param dims Numeric vector of dimensions for dimensionality reduction. Default: 10
#' @param clustmethod Character vector of clustering methods Default: clust.seurat
#' @param resolution Character vector of clustering resolutions. Default: 0.1, 1
#' 
#' @return List of factors containing clustering results from each pipeline run.
#' @export
#' 
#' @details
#' Valid filtering methods:
#'filt.default, filt.lenient, filt.mad, filt.pca, filt.pca2, filt.stringent, none
#'
#' Valid normalization methods:
#' norm.none, norm.none.scaled, norm.scnorm, norm.scnorm.scaled, norm.scran, norm.scran.scaled,
#' norm.sctransform, norm.scVI, norm.seurat, norm.seuratvst
#' 
#' Valid feature selection methods:
#' sel.deviance, sel.expr, sel.fromField, sel.vst
#' 
#' Valid dimensionality reduction method:
#' seurat.pca, serat.pca.noweight
#' 
#' @examples 
#' # Example 1: Run with default parameters
#' data(embryo)
#' res <- clustREval::runPipelineCombs(sce=embryo, outputPrefix="embryo")
#' res
#' 
#' # Example 2: Run with user-defined parameters
#' \dontrun{
#' data(embryo)
#' filtMethods <- c("filt.lenient", "filt.stringent")
#' res <- runPipelineCombs(sce=embryo, outputPrefix="embryo", filt=filtMethods)
#' res
#' }

runPipelineCombs <- function(sce, outputPrefix = "sce", doubletmethod =c("none"), filt=c("filt.default"), norm=c("norm.seurat"), sel=c("sel.vst"), selnb=2000, dr=c("seurat.pca"), clustmethod=c("clust.seurat"), dims=c(10), resolution=c(0.1, 1)){
  if(class(sce)[[1]] != "SingleCellExperiment"){
    stop("sce must be a SingleCellExperiment object.")
  }
  
  validFilt <- c("filt.default", "filt.lenient", "filt.mad", "filt.pca", "filt.pca2", "filt.stringent", "none")
  validNorm <- c("norm.none", "norm.none.scaled", "norm.scnorm", "norm.scnorm.scaled", "norm.scran", "norm.scran.scaled",
                 "norm.sctransform", "norm.scVI", "norm.seurat", "norm.seuratvst")
  validSel <- c("sel.deviance", "sel.expr", "sel.fromField", "sel.vst")
  drValid <- c("seurat.pca", "seurat.pca.noweight")
  
  if(!all(filt %in% validFilt)){
    stop("Filtering method specified not valid. Please specify a valid method.")
  }
  if(!all(norm %in% validNorm)){
    stop("Normalization method specified not valid. Please specify a valid method.")
  }
  if(!all(sel %in% validSel)){
    stop("Feature selection method specified not valid. Please specify a valid method.")
  }
  if(!all(dr %in% drValid)){
    stop("Dimensionality reduction method specified not valid. Please specify a valid method.")
  }
  
  if(any(resolution <= 0)){
    stop("Clustering resolutions must be greater than 0.")
  }
  if(any(dims <= 0)){
    stop("Dimensionaity reduction must have a dimension parameter greater than 0.")
  }
  
  alternatives <- list(
    doubletmethod=doubletmethod,
    filt=filt,
    norm=norm,
    sel=sel,
    selnb=selnb, 
    dr=dr,
    clustmethod=clustmethod,
    dims=dims,
    resolution=resolution
  )
  source(system.file("extdata", "scrna_alternatives.R", package="pipeComp"))
  
  sce <- scuttle::addPerCellQC(sce)
  sce <- scuttle::addPerFeatureQC(sce)
  # Throws a note since it's sourced from an external file, maybe copy/paste the function into this file?
  sce <- add_meta(sce)

  pct_counts_in_top_50_features <- sce$pct_counts_top_50_features
  sce$pct_counts_in_top_50_features <- pct_counts_in_top_50_features
  sce$pct_counts_top_50_features <- NULL
  sce$pct_counts_Mt <- sce$pct_Mt
  sce$pct_Mt <- NULL
  sce$pct_counts_in_top_20_features <- sce$percent.top_20
  sce$percent.top_20 <- NULL
  
  sce$phenoid <- rep(0, ncol(sce))

  sce <- c(sce=sce)
  pip_def <- pipeComp::scrna_pipeline(pipeClass = "seurat")
  evalDummy <- function(...){
    return(1)
  }
  pip_def@evaluation$doublet <- evalDummy
  pip_def@evaluation$filtering <- evalDummy
  pip_def@evaluation$normalization <- evalDummy
  pip_def@evaluation$selection <- evalDummy
  pip_def@evaluation$dimreduction <- evalDummy
  pip_def@evaluation$clustering <- evalDummy
  
  
  pip_def@aggregation$doublet <- evalDummy
  pip_def@aggregation$filtering <- evalDummy
  pip_def@aggregation$normalization <- evalDummy
  pip_def@aggregation$selection <- evalDummy
  pip_def@aggregation$dimreduction <- evalDummy
  pip_def@aggregation$clustering <- evalDummy
  
  outputPrefix <- paste("data", outputPrefix, sep="/")
  suppressWarnings(res <- pipeComp::runPipeline(sce, output.prefix = outputPrefix, alternatives, pip_def, nthreads=1, debug=FALSE))
  clustersPath <- paste0(outputPrefix, "res.sce.endOutputs.rds")
  
  print(clustersPath)
  clustRes <- readRDS(clustersPath)
  clustRes <- lapply(clustRes, labelled::remove_attributes, attributes="true.labels")
  return(clustRes)
}