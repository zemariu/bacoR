test_that("original BACO interpretation follows Diaz-Quijano categories", {
  a <- .baco_original_interpretation(1)
  b <- .baco_original_interpretation(1.5)
  c <- .baco_original_interpretation(0.5)
  d <- .baco_original_interpretation(-0.5)

  expect_identical(a$classification, "no_bias")
  expect_identical(b$classification, "overestimation")
  expect_identical(c$classification, "underestimation")
  expect_identical(d$classification, "inversion")
})

test_that("baco_effects stores interpretation when requested", {
  x <- baco_effects(
    effect_composite = 0.80,
    effect_critical = 0.60,
    se_log_composite = 0.08,
    se_log_critical = 0.12,
    correlation_log = 0.60,
    effect_measure = "RR",
    interpret = TRUE
  )

  expect_true(x$interpret_enabled)
  expect_identical(x$interpretation$baco$classification, "underestimation")
  expect_identical(
    x$interpretation$testing$conclusion,
    "evidence_of_discrepancy"
  )
})
