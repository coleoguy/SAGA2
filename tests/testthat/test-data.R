test_that("BanInf dataset loads and has correct structure", {
  data(BanInf)
  expect_s3_class(BanInf, "data.frame")
  expect_equal(ncol(BanInf), 7)
  expect_named(BanInf, c("cross", "mean", "SE", "sex", "enviro", "sire", "dam"))
  expect_true(nrow(BanInf) > 0)
})

test_that("BanInf has required parental lines", {
  data(BanInf)
  expect_true("P1" %in% BanInf$cross)
  expect_true("P2" %in% BanInf$cross)
})
