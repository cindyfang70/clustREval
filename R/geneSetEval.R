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
#' @param gmtPathway File path pointing to a .gmt file for gene set analysis.
#'
#' @return Tibble containing gene set enrichment results for all clusters.
#' @export
#'
#' @examples
#' # Example:
#' \dontrun{
#' data(embryo)
#' data(embryoClusts)
#' # Compute pathways for the first clustering output by accessing raw data provided for pathways file
#' GMTPath <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "clustREval")
#' geneSetRes <- geneSetEval(embryo, embryoClusts[[1]], GMTPath)
#' geneSetRes
#' }

geneSetEval <- function(sce, clusters, gmtPathway){
  if (!file.exists(gmtPathway)){
    stop("File not found. Please provide a valid filepath for the hallmark pathways.")
  }
  if(!grepl("*\\.gmt$",gmtPathway)){
    stop("File format not valid, please provide a .gmt file.")
  }
  if(class(sce)[[1]] != "SingleCellExperiment"){
    stop("sce must be a SingleCellExperiment object.")
  }
  if(class(clusters)[[1]] != "factor"){
    stop("Clustering results must be provided as a factor.")
  }
  if(length(clusters) > dim(sce)[[2]]){
    stop("There are more cells in the clustering output than in the SingleCellExperiment. Please ensure that all cells
         in the clustering output are from the original experiment.")
  }
  if(!any(rownames(as.matrix(clusters)) %in% colnames(sce))){
    stop("Clustered cells must be a subset of cells in the SingleCellExperiment object.")
  }

  pathways.hallmark <- fgsea::gmtPathways(gmtPathway)
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
    res <- dplyr::inner_join(tibble::as_tibble(logfc, rownames="row"), ens2symbol, by=c("row"="ENSEMBL"))
    res2 <- res %>%
      dplyr::select("SYMBOL", "value") %>%
      stats::na.omit() %>%
      dplyr::distinct() %>%
      dplyr::group_by("SYMBOL")
    ranks <- tibble::deframe(res2)
    ranks
  })
  
  gseas <- lapply(mapped_logFCs, function(mapped_logfc){
    fgseaRes <- fgsea::fgsea(mapped_logfc, pathways = pathways.hallmark)
    fgseaResTidy <- fgseaRes %>%
      tibble::as_tibble() %>%
      dplyr::arrange(dplyr::desc(NES))
    #mean(abs(fgseaResTidy$ES))
  })
  return(gseas)
}