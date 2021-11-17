# test_that("Results are of expected length", {
#   res1 <- runPipelineCombs(embryo)
#   expect_equal(length(res1), 2)
#   res2 <- runPipelineCombs(embryo, filt=c("filt.stringent", "filt.default"))
#   expect_equal(length(res2),4)
# })
# 
# test_that("Handles invalid input", {
#   expect_error(runPipelineCombs(as.matrix(embryo)))
#   expect_error(runPipelineCombs(embryo, filt=c("asdf")))
#   expect_error(runPipelineCombs(embryo, dims=c(-1)))
# })