library(pipeComp)


#' runPipelineCombs
#'
#'Run user-defined clustering pipeline combinations on a user-specificed scRNA-seq dataset.
#' @param sce SingleCellExperiment object with counts matrix and phenoid colData entry.
#' @param filt list of filtering methods
#' @param norm list of normalization methods
#' @param resolution list of clustering resolutions
#'
#' @return end outputs object
#' @export

runPipelineCombs <- function(sce, filt, norm, resolution){
  
  alternatives <- list(
    filt=filt,
    norm=norm,
    resolution=resolution
  )
  source(system.file("extdata", "scrna_alternatives.R", package="pipeComp"))
  pip_def <- pipeComp::scrna_pipeline(pipeClass = "seurat")
  res <- pipeComp::runPipeline(sce, alternatives, pip_def, nthreads=1)
}