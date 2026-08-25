# The exhaustive model-scoring pass uses lm.wfit for speed (saga2-12).
# This test keeps a small glm-based reference implementation and asserts
# that AICc values and weights from LCA match the glm path exactly on the
# shipped BanInf data.

test_that("lm.wfit scoring pass matches a glm reference implementation", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)

  # rebuild the model space exactly as AnalyzeModels does
  cmat <- result$cmatrix[, -1]
  mm <- as.matrix(cmat)
  storage.mode(mm) <- "numeric"
  pos.cols <- 2:ncol(mm)
  max.par <- min(nrow(mm) - 2, length(pos.cols))
  eqns <- list()
  counter <- 1
  for (i in 1:max.par) {
    foo <- combn(pos.cols, i)
    for (j in 1:ncol(foo)) {
      eqns[[counter]] <- as.vector(foo[, j])
      counter <- counter + 1
    }
  }

  # glm reference: fit every model, drop singular fits, compute AICc
  weights <- BanInf$SE^-2
  n <- nrow(BanInf)
  ref.aic <- rep(NA_real_, length(eqns))
  ref.rank <- rep(NA_real_, length(eqns))
  for (i in seq_along(eqns)) {
    g <- glm(BanInf$mean ~ mm[, c(1, eqns[[i]])], weights = weights)
    if (sum(is.na(g$coef)) < 2) {
      ref.aic[i] <- g$aic
      ref.rank[i] <- g$rank
    }
  }
  kept <- !is.na(ref.aic)
  K <- ref.rank[kept] + 1
  ref.aicc <- ref.aic[kept] + (2 * K * (K + 1)) / (n - K - 1)
  ref.daicc <- ref.aicc - min(ref.aicc)

  # same number of retained models and identical AICc-based results
  expect_equal(length(result$all.models), sum(kept))
  expect_equal(as.numeric(result$daicc), as.numeric(ref.daicc))

  # Akaike weights identical too
  ref.w <- exp(-0.5 * ref.daicc) / sum(exp(-0.5 * ref.daicc))
  w <- exp(-0.5 * result$daicc) / sum(exp(-0.5 * result$daicc))
  expect_equal(as.numeric(w), as.numeric(ref.w))
})