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
  gseas$pathway <- gsub("HALLMARK_", "", gseas$pathway)
  p <- ggplot2::ggplot(gseas, ggplot2::aes(reorder(pathway, NES), NES)) +
      ggplot2::geom_col(ggplot2::aes(fill=padj<0.05)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x="Pathway", y="Normalized Enrichment Score",
        title=sprintf("Hallmark pathways NES from GSEA for Cluster"))+
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=3))
  return(p)

}