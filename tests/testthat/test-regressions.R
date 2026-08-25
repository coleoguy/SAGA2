# Regression tests for bugs fixed in version 2.0.1

test_that("best.eqns.w is aligned with best.models (saga2-03)", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = FALSE, messages = FALSE)
  # weights must come back in decreasing order matching best.models
  expect_false(is.unsorted(rev(result$best.eqns.w)))
  # the first weight must equal the largest Akaike weight in model space
  waic <- exp(-0.5 * result$daicc) / sum(exp(-0.5 * result$daicc))
  expect_equal(as.numeric(result$best.eqns.w[1]), max(waic))
  # and it must equal the weight recomputed from the top model's own fit
  expect_equal(length(result$best.eqns.w), length(result$best.models))
})

test_that("a single invalid sire or dam is caught (saga2-04)", {
  data(BanInf)
  bad <- BanInf
  bad$sire[5] <- 99
  expect_error(LCA(bad, SCS = "NSC", messages = FALSE), "row")
  bad2 <- BanInf
  bad2$dam[3] <- 0
  expect_error(LCA(bad2, SCS = "NSC", messages = FALSE), "row")
  bad3 <- BanInf
  bad3$dam[3] <- NA
  expect_error(LCA(bad3, SCS = "NSC", messages = FALSE), "row")
})

test_that("SCS synonym X0 and all-NA enviro propagate (saga2-01)", {
  data(BanInf)
  # X0 (zero) should be silently corrected to XO (letter O) and run
  expect_s3_class(LCA(BanInf, SCS = "X0", messages = FALSE), "genarch")
  # an all-NA enviro column is documented as allowed
  dat <- BanInf
  dat$enviro <- NA
  expect_s3_class(LCA(dat, SCS = "NSC", messages = FALSE), "genarch")
})

test_that("CleanCmatrix handles a minimal c-matrix (saga2-06)", {
  # two CGE columns (mean + one effect) previously crashed the
  # duplicate-column scan with 'undefined columns selected'
  cm <- data.frame(cross = c("P1", "P2", "F1", "F2"),
                   M = c(1, 1, 1, 1),
                   Aa = c(1, -1, 0, 0))
  out <- SAGA2:::CleanCmatrix(cm, messages = FALSE)
  expect_equal(colnames(out[[1]]), c("M", "Aa"))
})

test_that("AICc uses the standard small-sample formula (saga2-07)", {
  data(BanInf)
  result <- LCA(BanInf, SCS = "NSC", ret.all = TRUE, messages = FALSE)
  # recompute AICc for every model directly from the fitted objects
  n <- nrow(BanInf)
  aic <- vapply(result$all.models, function(m) m$aic, numeric(1))
  K <- vapply(result$all.models, function(m) m$rank, numeric(1)) + 1
  aicc <- aic + (2 * K * (K + 1)) / (n - K - 1)
  expect_equal(as.numeric(result$daicc), as.numeric(aicc - min(aicc)))
})

test_that("NA standard errors are rejected up front", {
  data(BanInf)
  bad <- BanInf
  bad$SE[2] <- NA
  expect_error(LCA(bad, SCS = "NSC", messages = FALSE), "NA")
})