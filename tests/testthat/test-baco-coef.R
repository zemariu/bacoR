test_that("baco_coef reproduces the worked adjusted example", {
  x <- baco_coef(
    beta_composite = log(0.80),
    beta_critical = log(0.60),
    se_composite = 0.08,
    se_critical = 0.12,
    covariance = 0.00576,
    effect_measure = "RR"
  )

  expect_equal(x$baco$estimate, log(0.80) / log(0.60), tolerance = 1e-12)
  expect_equal(x$discrepancy$delta, log(0.80) - log(0.60), tolerance = 1e-12)
  expect_equal(x$discrepancy$effect_ratio, 0.80 / 0.60, tolerance = 1e-12)
  expect_equal(x$joint$correlation, 0.60, tolerance = 1e-12)
  expect_identical(x$baco$fieller$set_type, "bounded")
  expect_equal(x$baco$fieller$intervals$lower, 0.1692776, tolerance = 1e-5)
  expect_equal(x$baco$fieller$intervals$upper, 0.7241962, tolerance = 1e-5)
})

test_that("baco_coef rejects impossible covariance", {
  expect_error(
    baco_coef(
      beta_composite = -0.2,
      beta_critical = -0.5,
      se_composite = 0.08,
      se_critical = 0.12,
      covariance = 0.02
    ),
    "incompatible"
  )
})

test_that("baco_coef identifies a weak denominator", {
  x <- baco_coef(
    beta_composite = -0.10,
    beta_critical = -0.01,
    se_composite = 0.04,
    se_critical = 0.05,
    covariance = 0.001,
    effect_measure = "RR"
  )

  expect_true(x$diagnostics$weak_denominator)
})
