test_that("Fieller returns the known bounded example", {
  out <- baco_fieller(
    beta_composite = log(0.80),
    beta_critical = log(0.60),
    var_composite = 0.08^2,
    var_critical = 0.12^2,
    covariance = 0.00576
  )

  expect_equal(out$estimate, log(0.80) / log(0.60), tolerance = 1e-12)
  expect_identical(out$set_type, "bounded")
  expect_equal(out$intervals$lower, 0.1692776, tolerance = 1e-5)
  expect_equal(out$intervals$upper, 0.7241962, tolerance = 1e-5)
  expect_false(out$weak_denominator)
})

test_that("Fieller identifies weak denominator settings", {
  out <- baco_fieller(
    beta_composite = log(0.98),
    beta_critical = log(1.01),
    var_composite = 0.02^2,
    var_critical = 0.05^2,
    covariance = 0.0005
  )
  expect_true(out$weak_denominator)
  expect_true(out$set_type %in% c("disjoint_unbounded", "all_real", "half_line"))
})
