library(cluster)
library(Seurat)
library(clValid)

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
  var.names <- head(Seurat::VariableFeatures(var.genes), 500)
  filt_sce <- filt_sce[var.names,]
  
  clusters <- as.numeric(clusters)
  dunnIndex <- clValid::dunn(clusters=clusters, Data=as.matrix(counts(filt_sce)))
  dissim <- cluster::daisy(t(as.matrix(counts(filt_sce))))
  silhouetteIndex <- cluster::silhouette(as.matrix(clusters), dist=dissim)

  metrics <- cbind(dunnIndex, mean(silhouetteIndex[,3]))
  return(metrics)
}