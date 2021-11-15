library(Seurat)
library(gridExtra)

#' Title
#'
#' @param sce A single cell experiment object
#' @param clusters Clustering results for the SingleCellExperiment
#' @param gseas A tibble containing enrichment values for each cluster computed from geneSetEval
#'
#' @return a plot
#' @export
#'
#' @examples
plotGeneSetEval <- function(sce, clusters, gseas){
  plots <- c()
  for(i in 1:length(gseas)){
    gseas[[i]]$pathway <- gsub("HALLMARK_", "", gseas[[i]]$pathway)
    p <- ggplot(gseas[[i]], aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj<0.05)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
        title=sprintf("Hallmark pathways NES from GSEA for Cluster %s", names(gseas)[[i]]))+
      theme(axis.text.y=element_text(size=3))
    plots <- c(plots, list(p))
  }
  n <- length(plots)
  nCol <- floor(sqrt(n))
  do.call("grid.arrange", c(plots=plots, ncol=nCol))
}