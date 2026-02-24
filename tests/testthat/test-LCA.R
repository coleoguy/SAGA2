test_that("LCA returns a genarch object", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
  expect_s3_class(result, "genarch")
})

test_that("LCA result has expected components", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
  expect_named(result, c("all.models", "best.models", "best.eqns.w",
                          "estimates", "daicc", "varimp", "cmatrix"))
})

test_that("LCA estimates matrix has correct dimensions", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
  expect_equal(nrow(result$estimates), 2)
  expect_true(ncol(result$estimates) > 1)
  expect_equal(rownames(result$estimates),
               c("Model Weighted Average", "Unconditional Standard Error"))
})

test_that("LCA with ret.all returns model objects", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)
  expect_true(!is.null(result$all.models))
  expect_true(length(result$all.models) > 0)
})

test_that("LCA without ret.all returns NULL for all.models", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = FALSE, messages = FALSE)
  expect_null(result$all.models)
})

test_that("LCA varimp scores are between 0 and 1", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
  vi_scores <- as.numeric(result$varimp[, 2])
  expect_true(all(vi_scores >= 0 & vi_scores <= 1))
})

test_that("LCA drop.pars removes specified CGEs", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", drop.pars = c("Ca"),
                messages = FALSE)
  expect_false("Ca" %in% colnames(result$estimates))
})

test_that("LCA max.pars limits model complexity", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", max.pars = 2, messages = FALSE)
  expect_s3_class(result, "genarch")
})
