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
#' @param doubletmethod Character list of doublet detection methods. Default: none
#' @param filt Character list of filtering methods. Default: filt.default
#' @param norm Character list of normalization methods Default: norm.seurat
#' @param sel Character list of selection methods. Default: sel.vst
#' @param selnb Numeric list of number of genes to select. Default: 2000
#' @param dr Character list of dimensionality reduction methods. Default: seurat.pca
#' @param dims Numeric list of dimensions for dimensionality reduction. Default: 10
#' @param clustmethod Character list of clustering methods Default: clust.seurat
#' @param resolution Character list of clustering resolutions. Default: 0.1, 0.2
#' 
#' @return List of factors containing clustering results from each pipeline run.
#' @export
#' 
#' @examples 
#' # Example 1: Run with default parameters
#' data(embryo)
#' res <- runPipelineCombs(sce=embryo, outputPrefix="results/embryo")
#' res
#' 
#' # Example 2: Run with user-defined parameters
#' data(embryo)
#' filtMethods <- c("filt.lenient", "filt.stringent")
#' res <- runPipelineCombs(sce=embryo, outputPrefix="results/embryo", filt=filtMethods)
#' res

runPipelineCombs <- function(sce, outputPrefix = "sce", doubletmethod =c("none"), filt=c("filt.default"), norm=c("norm.seurat"), sel=c("sel.vst"), selnb=2000, dr=c("seurat.pca"), clustmethod=c("clust.seurat"), dims=c(10), resolution=c(0.1, 0.2)){
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
  res <- pipeComp::runPipeline(sce, output.prefix = outputPrefix, alternatives, pip_def, nthreads=1, debug=TRUE)
  clustersPath <- paste0(outputPrefix, "res.sce.endOutputs.rds")
  
  print(clustersPath)
  clustRes <- readRDS(clustersPath)
  clustRes <- lapply(clustRes, labelled::remove_attributes, attributes="true.labels")
  return(clustRes)
}