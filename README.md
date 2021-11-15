
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ClustREval

<!-- badges: start -->

<!-- badges: end -->

## Description

ClustREval is a package for evaluating the performance of different
clustering pipelines on scRNA-seq data using unsupervised metrics and
gene set enrichment analysis. Clustering is a powerful tool for helping
researchers to detect cellular heterogeniety. However, clustering
performance is highly dependent on parameters used in the clustering
pipeline, for which there are no systematic recommendations. This
package allows users to compute clustering results from various
clustering pipelines defined by user-specified parameters. Clustering
results can then be evaluated using comparing unsupervised clustering
metrics and differential gene expression between results.

## Installation

To install the latest version of the package:

``` r
require("devtools")
devtools::install_github("cindyfang70/clustREval", build_vignettes = TRUE)
library("clustREval")
```

## Overview

``` r
ls("package:clustREval")
data(package = "clustREval") # optional
```

``` r
browseVignettes("clustREval")
```

The package tree structure is provided below (optional).

## Contributions

The contributor for this package is Xin Zhi Fang.

## References

1.  Germain, P. L., Sonrel, A., & Robinson, M. D. (2020). pipeComp, a
    general framework for the evaluation of computational pipelines,
    reveals performant single cell RNA-seq preprocessing tools. Genome
    Biology, 21(1). <https://doi.org/10.1186/s13059-020-02136-7>

2.  Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of
    fold change and dispersion for RNA-seq data with DESeq2. Genome
    Biology, 15(12). <https://doi.org/10.1186/s13059-014-0550-8>

3.  Amezquita, R. A., Lun, A. T. L., Becht, E., Carey, V. J., Carpp, L.
    N., Geistlinger, L., Marini, F., Rue-Albrecht, K., Risso, D.,
    Soneson, C., Waldron, L., Pagès, H., Smith, M. L., Huber, W.,
    Morgan, M., Gottardo, R., & Hicks, S. C. (2019). Orchestrating
    single-cell analysis with Bioconductor. Nature Methods, 17(2),
    137–145. <https://doi.org/10.1038/s41592-019-0654-x>

4.  Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N.,
    & Sergushichev, A. (2016). Fast gene set enrichment analysis.
    BioRxiv. Published. <https://doi.org/10.1101/060012>

## Acknowledgements

This package was developed for BCB410H: Applied Bioinformatics,
University of Toronto, Toronto, CANADA, 2019-2021.
