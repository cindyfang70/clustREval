test_that("Handles invalid input", {
  data(embryo)
  data(embryoClusts)
  expect_error(computeUnsupervisedMetrics(embryo, as.matrix(embryoClusts[[1]])))
  
  fakeClusters <- as.factor(rep(1, 125))
  expect_error(computeUnsupervisedMetrics(embryo, fakeClusters))
  
  clusters2 <- as.factor(c(embryoClusts[[1]], asdf=1))
  expect_error(computeUnsupervisedMetrics(embryo, clusters2))
  
})

test_that("Results are of expected length", {
  data(embryo)
  data(embryoClusts)
  metrics <- computeUnsupervisedMetrics(embryo, embryoClusts[[1]])
  expect_equal(length(metrics),2)
})

