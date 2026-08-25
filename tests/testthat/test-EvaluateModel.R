test_that("EvaluateModel returns expected structure", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)
  eval_result <- EvaluateModel(result, model = 1)
  expect_type(eval_result, "list")
  expect_named(eval_result, c("Estimates", "Estiamtes", "P.val",
                              "Bonf.Corr.Alpha"))
  # the deprecated alias must match the new element
  expect_identical(eval_result$Estimates, eval_result$Estiamtes)
})

test_that("EvaluateModel Bonferroni alpha uses the number of models", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)
  eval_result <- EvaluateModel(result, model = 1)
  expect_equal(eval_result$Bonf.Corr.Alpha, 0.05 / length(result$all.models))
})

test_that("EvaluateModel errors informatively without all.models", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = FALSE, messages = FALSE)
  expect_error(EvaluateModel(result, model = 1), "ret.all")
})

test_that("EvaluateModel returns numeric estimates", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)
  eval_result <- EvaluateModel(result, model = 1)
  expect_true(is.numeric(eval_result$Estiamtes))
  expect_true(is.numeric(eval_result$P.val))
  expect_true(is.numeric(eval_result$Bonf.Corr.Alpha))
})
