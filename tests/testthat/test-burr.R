test_that("Burr3 CDF and quantile invert", {
  p <- c(0.01, 0.2, 0.5, 0.9, 0.99)
  x <- qburr3(p, c=2, k=3, s=1)
  expect_equal(pburr3(x, c=2, k=3, s=1), p, tolerance = 1e-8)
})
