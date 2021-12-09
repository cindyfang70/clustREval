#' Tracing pluripotency of human early embryos and embryonic stem cells by single cell RNA-seq
#' 
#'
#' Gene expression of 124 early embryo and embryonic stem cells from scRNA-seq. Sourced from the
#' EMBL-EBI Single Cell Expression Atlas.
#' 
#' Yan L, Yang M, Guo H, et al. Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nature Structural & Molecular Biology. 2013 Sep;20(9):1131-1139. DOI: 10.1038/nsmb.2660. PMID: 23934149.
#' @format A data frame with 26073 rows and 124 columns:
#' @source \url{https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-36552/results/tsne}
"embryo"

#' Embryo data clustering results from two pipelines.
#' 
#'
#' Clustering results of the embryo scRNA-seq data from two pipelines. Pipelines computed are the default pipeline settings from runPipelineCombs.
#' Structured as a list with two entries, each entry is a factor that indicates the cluster that each of the cells in the embryo datasets belong to.
#' @format A list with two levels
"embryoClusts"

#' GSEA results for embryo data clustered with first pipeline option.
#' 
#'
# GSEA results for the four clusters outputted by the pipeline.
#' @format A list with 4 levels.
"gsea"


