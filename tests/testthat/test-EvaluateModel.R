test_that("EvaluateModel returns expected structure", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)
  eval_result <- EvaluateModel(result, model = 1)
  expect_type(eval_result, "list")
  expect_named(eval_result, c("Estiamtes", "P.val", "Bonf.Corr.Alpha"))
})

test_that("EvaluateModel returns numeric estimates", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)
  eval_result <- EvaluateModel(result, model = 1)
  expect_true(is.numeric(eval_result$Estiamtes))
  expect_true(is.numeric(eval_result$P.val))
  expect_true(is.numeric(eval_result$Bonf.Corr.Alpha))
})
