library(cluster)
library(Seurat)
library(clValid)
library(SingleCellExperiment)

#' computeUnsupervisedMetrics
#'
#' A function that computes Dunn Index and and mean Silhouette score for a clustering output
#'
#' @param sce A SingleCellExperiment object with raw counts in the counts slot.
#' @param clusters Factor representing clustering results for the SingleCellExperiment object.
#'
#' @return A list with the Dunn Index and mean Silhouette score as elements
#' @export
#'
#' @examples
#' # Example:
#' \dontrun{
#' data(embryo)
#' data(embryoClusts)
#' # Compute metrics for the first clustering output on the embryo dataset
#' metrics <- computeUnsupervisedMetrics(embryo, embryoClusts[[1]])
#' }

#' 
computeUnsupervisedMetrics <- function(sce, clusters){
  
  if(class(sce)[[1]] != "SingleCellExperiment"){
    stop("sce must be a SingleCellExperiment object.")
  }
  if(class(clusters)[[1]] != "factor"){
    stop("Clustering results must be provided as a factor.")
  }
  # if(length(clusters) > dim(sce)[[2]]){
  #   stop("There are more cells in the clustering output than in the SingleCellExperiment. Please ensure that all cells
  #        in the clustering output are from the original experiment.")
  # }
  if(!any(rownames(as.matrix(clusters)) %in% colnames(sce))){
    stop("Clustered cells must be a subset of cells in the SingleCellExperiment object.")
  }
  filt_sce <- sce[,rownames(as.matrix(clusters))]
  filt_sce <- scuttle::logNormCounts(filt_sce)


  seu <- Seurat::as.Seurat(filt_sce, counts="counts", data="logcounts")
  var.genes <- Seurat::FindVariableFeatures(seu)
  var.names <- utils::head(Seurat::VariableFeatures(var.genes), 500)
  filt_sce <- filt_sce[var.names,]
  
  clusters <- as.numeric(clusters)
  countsMat <- as.matrix(SingleCellExperiment::counts(filt_sce))
  dunnIndex <- clValid::dunn(clusters=clusters, Data=countsMat)
  dissim <- cluster::daisy(t(countsMat))
  silhouetteIndex <- cluster::silhouette(as.matrix(clusters), dist=dissim)

  metrics <- c(dunnIndex=dunnIndex, meanSilhouette=mean(silhouetteIndex[,3]))
  return(metrics)
}