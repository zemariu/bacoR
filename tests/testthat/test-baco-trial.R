test_that("baco_trial is structurally valid and nested", {
  expect_s3_class(baco_trial, "data.frame")
  expect_equal(nrow(baco_trial), 2500)
  expect_true(all(c(
    "id", "treatment", "age", "sex", "comorbidity",
    "composite", "critical"
  ) %in% names(baco_trial)))
  expect_true(all(baco_trial$critical <= baco_trial$composite))

  pars <- attr(baco_trial, "data_generating_parameters")
  expect_equal(pars$rr_composite, 0.80)
  expect_equal(pars$rr_critical, 0.60)
  expect_equal(pars$baco, log(0.80) / log(0.60), tolerance = 1e-12)
})
