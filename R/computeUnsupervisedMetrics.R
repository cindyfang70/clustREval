library(cluster)
library(Seurat)
library(clValid)
library(SingleCellExperiment)

#' computeUnsupervisedMetrics
#'
#' @param sce A SingleCellExperiment object
#' @param clusters Clustering results for the SingleCellExperiment object
#'
#' @return unsupervised metrics
#' @export
#'
#' @examples
computeUnsupervisedMetrics <- function(sce, clusters){
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

  metrics <- cbind(dunnIndex, mean(silhouetteIndex[,3]))
  return(metrics)
}