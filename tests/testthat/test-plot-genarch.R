test_that("plot.genarch runs without error", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
  expect_no_error(plot(result, min.vi = 0.3))
})
