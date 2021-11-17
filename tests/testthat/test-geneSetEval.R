test_that("Handles invalid input",{
  GMTPath <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "clustREval")
  expect_error(geneSetEval(embryo, embryoClusts[[1]], "asdf.gmt"))
  expect_error(geneSetEval(embryo, embryoClusts[[1]], "README.md"))
  expect_error(geneSetEval(embryo, as.matrix(embryoClusts[[1]]), GMTPath))
  
  fakeClusters <- as.factor(rep(1, dim(embryo)[[2]]+1))
  expect_error(geneSetEval(embryo, fakeClusters, GMTPath))
  
  clusters2 <- as.factor(c(embryoClusts[[1]], asdf=1))
  expect_error(geneSetEval(embryo, clusters2, GMTPath))
})

test_that("Results are of expected length", {
  GMTPath <- system.file("extdata", "h.all.v7.4.symbols.gmt", package = "clustREval")
  gsea <- geneSetEval(embryo, embryoClusts[[1]], GMTPath)
  nclusts <- length(unique(embryoClusts[[1]]))
  expect_equal(nclusts, length(gsea))
})
