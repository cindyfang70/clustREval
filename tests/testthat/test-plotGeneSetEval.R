test_that("Handles invalid input", {
  testGsea <- c("pathway")
  expect_error(plotGeneSetEval(testGsea))
})

