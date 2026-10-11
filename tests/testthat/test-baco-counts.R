test_that("regular nested counts reproduce expected BACO and Fisher result", {
  out <- baco_counts(
    n_treatment = 1000,
    n_control = 1000,
    composite_treatment = 150,
    composite_control = 200,
    critical_treatment = 30,
    critical_control = 60
  )

  expect_equal(out$effects$rr_composite, 0.75)
  expect_equal(out$effects$rr_critical, 0.50)
  expect_equal(out$baco$estimate, log(0.75) / log(0.50), tolerance = 1e-12)
  expect_equal(out$discrepancy$delta, log(0.75) - log(0.50), tolerance = 1e-12)
  expect_equal(out$discrepancy$qrr, (30 / 150) / (60 / 200), tolerance = 1e-12)
  expect_lt(out$exact_test$fisher_p, 0.05)
})

test_that("one zero critical cell leaves Fisher available but BACO non-estimable", {
  out <- baco_counts(
    n_treatment = 500,
    n_control = 500,
    composite_treatment = 60,
    composite_control = 80,
    critical_treatment = 0,
    critical_control = 4
  )

  expect_false(out$diagnostics$baco_estimable)
  expect_true(is.finite(out$exact_test$fisher_p))
  expect_true(out$diagnostics$zero_critical_any_arm)
})

test_that("zero critical events in both groups is reported as no information", {
  out <- baco_counts(
    n_treatment = 500,
    n_control = 500,
    composite_treatment = 50,
    composite_control = 55,
    critical_treatment = 0,
    critical_control = 0
  )

  expect_false(out$diagnostics$baco_estimable)
  expect_equal(out$exact_test$fisher_p, 1)
  expect_true(out$diagnostics$zero_critical_both_arms)
  expect_match(out$exact_test$note, "absence of critical-event information")
})

test_that("invalid non-nested counts are rejected", {
  expect_error(
    baco_counts(100, 100, 10, 10, 11, 2),
    "nested"
  )
})
