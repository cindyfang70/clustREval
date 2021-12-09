#' Plot Enrichment Scores for Each Pathway in a Cluster
#'
#'  A function that plots gene set enrichment analysis results for a cluster in
#' order to visualize biological differences between clusters.
#'
#' @param gseas A list containing enrichment values for each pathway in a cluster
#' @param clustname A string or integer containing the name of the cluster. Used in the title of the plot.
#'
#' @return A plot for one cluster from a clustering output showing the enriched pathways for that cluster.
#' 
#' @export
#'
#' @examples
#' # Example:
#' \dontrun{
#' data(embryo)
#' data(embryoClusts)
#' # Compute pathways for the first clustering output
#' data(embryo)
#' data(embryoClusts)
#' GMTPath <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "clustREval")
#' geneSetRes <- geneSetEval(embryo, embryoClusts[[1]], GMTPath)
#' plotGeneSetEval(geneSetRes[[1]])
#' }
#' 
#' @import ggplot2

plotGeneSetEval <- function(gseas, clustname=""){
  
  if(length(gseas) < 1){
    stop("You must provide a tibble with at least one hallmark pathway.")
  }
  
  gseas$pathway <- gsub("HALLMARK_", "", gseas$pathway)
  p <- ggplot2::ggplot(gseas, ggplot2::aes(reorder(pathway, NES), NES)) +
      ggplot2::geom_col(ggplot2::aes(fill=padj<0.05)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x="Pathway", y="Normalized Enrichment Score",
        title=sprintf("Hallmark pathways NES from GSEA for Cluster %s", clustname))+
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=3))
  return(p)
}
#[END]