library(pipeComp)


#' runPipelineCombs
#'
#'Run user-defined clustering pipeline combinations on a user-specificed scRNA-seq dataset.
#'
#' @param sce SingleCellExperiment object with counts matrix and phenoid colData entry.
#' @param filt list of filtering methods
#' @param norm list of normalization methods
#' @param resolution list of clustering resolutions
#' @param doubletmethod list of doublet detection methods
#' @param sel list of selection methods
#' @param selnb list of number of genes to select
#' @param dr list of dimensionality reduction methods
#' @param clustmethod list of clustering methods
#' @param dims list of dimensions for dimensionality reduction
#'
#' @return end outputs object
#' @export

runPipelineCombs <- function(sce, doubletmethod =c("none"), filt=c("default"), norm=c("norm.seurat"), sel=c("sel.vst"), selnb=2000, dr=c("seurat.pca"), clustmethod=c("clust.seurat"), dims=c(10), resolution=c(0.1, 0.2)){
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
  pip_def <- pipeComp::scrna_pipeline(pipeClass = "seurat")
  res <- pipeComp::runPipeline(sce, alternatives, pip_def, nthreads=1)
}