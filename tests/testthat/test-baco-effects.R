test_that("baco_effects logs adjusted effects internally", {
  x <- baco_effects(
    effect_composite = 0.80,
    effect_critical = 0.60,
    se_log_composite = 0.08,
    se_log_critical = 0.12,
    correlation_log = 0.60,
    effect_measure = "RR"
  )

  expect_equal(x$baco$estimate, log(0.80) / log(0.60), tolerance = 1e-12)
  expect_equal(x$joint$covariance, 0.00576, tolerance = 1e-12)
  expect_equal(x$discrepancy$effect_ratio, 0.80 / 0.60, tolerance = 1e-12)
})
