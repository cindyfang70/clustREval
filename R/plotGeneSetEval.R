library(Seurat)
library(gridExtra)
library(ggplot2)

#' plotGeneSetEval
#'
#'  A function that plots gene set enrichment analysis results for each cluster in
#' order to visualize biological differences between clusters.
#'
#' @param gseas A list containing enrichment values for each cluster computed from geneSetEval
#'
#' @return A plot for each clustering output showing the enriched pathways for each.
#' @export
#'
#' @examples
#' # Example:
#' \dontrun{
#' data(embryo)
#' data(embryoClusts)
#' # Compute pathways for the first clustering output
#' geneSetRes <- geneSetEval(embryo, embryoClusts[[1]])
#' plotGeneSetEval(geneSetRes)
#' }

plotGeneSetEval <- function(gseas){
  plots <- c()
  suppressPackageStartupMessages(require(ggplot2))
  for(i in 1:length(gseas)){
    gseas[[i]]$pathway <- gsub("HALLMARK_", "", gseas[[i]]$pathway)
    p <- ggplot2::ggplot(gseas[[i]], aes(reorder(pathway, NES), NES)) +
      ggplot2::geom_col(aes(fill=padj<0.05)) +
      ggplot2::coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
        title=sprintf("Hallmark pathways NES from GSEA for Cluster %s", names(gseas)[[i]]))+
      theme(axis.text.y=element_text(size=3))
    plots <- c(plots, list(p))
  }
  n <- length(plots)
  nCol <- floor(sqrt(n))
  suppressPackageStartupMessages(require(gridExtra))
  do.call("grid.arrange", c(plots=plots, ncol=nCol))
}