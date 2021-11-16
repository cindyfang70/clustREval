library(scran)
library(scuttle)
library(fgsea) # BioConductor
library(org.Hs.eg.db) # BioConductor
library(AnnotationDbi) # BioConductor
library(tidyverse)
library(Seurat)
library(dplyr)
library(tibble)

#' geneSetEval
#' 
#' A function that uses MSigDB Hallmark Pathways to compute gene set enrichment analysis for each clustering output 
#' given in clusters.
#'
#' @param sce A SingleCellExperiment object with raw counts in the counts slot.
#' @param clusters Factor representing clustering results from a pipeline for the SingleCellExperiment.
#'
#' @return Tibble containing gene set enrichment results for all clusters.
#' @export
#'
#' @examples
#' # Example:
#' \dontrun{
#' data(embryo)
#' data(embryoClusts)
#' # Compute pathways for the first clustering output
#' geneSetRes <- geneSetEval(embryo, embryoClusts[[1]])
#' geneSetRes
#' }

geneSetEval <- function(sce, clusters){
  pathways.hallmark <- fgsea::gmtPathways("h.all.v7.4.symbols.gmt")
  filt_sce <- sce[,rownames(as.matrix(clusters))]
  filt_sce <- scuttle::logNormCounts(filt_sce)
  
  fms <- scran::findMarkers(filt_sce, groups = clusters)
  logFCs <- lapply(fms,
                   function(fm){
                     print(fm)
                     logfc <- fm[[5]]
                     names(logfc) <- rownames(fm)
                     logfc
                   })
  
  mapped_logFCs <- lapply(logFCs, function(logfc){
    ens2symbol <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                        key=names(logfc),
                                        columns="SYMBOL",
                                        keytype="ENSEMBL")
    ens2symbol <- tibble::as_tibble(ens2symbol)
    res <- dplyr::inner_join(as_tibble(logfc, rownames="row"), ens2symbol, by=c("row"="ENSEMBL"))
    res2 <- res %>%
      dplyr::select("SYMBOL", value) %>%
      na.omit() %>%
      distinct() %>%
      group_by("SYMBOL")
    ranks <- deframe(res2)
    ranks
  })
  
  gseas <- lapply(mapped_logFCs, function(mapped_logfc){
    fgseaRes <- fgsea(mapped_logfc, pathways = pathways.hallmark)
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    #mean(abs(fgseaResTidy$ES))
  })
  return(gseas)
}