test_that("VisModelSpace runs without error", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
  expect_no_error(VisModelSpace(result))
})
