test_that("plotObserved runs without error", {
  data(BanInf)
  expect_no_error(plotObserved(BanInf))
})

test_that("plotObserved runs without SE bars", {
  data(BanInf)
  expect_no_error(plotObserved(BanInf, SE = FALSE))
})
