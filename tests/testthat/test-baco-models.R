test_that("baco_models extracts fitted GLM effects and joint covariance", {
  set.seed(101)
  n <- 1500
  treatment <- rbinom(n, 1, 0.5)
  age <- rnorm(n, 50, 10)
  comorbidity <- rbinom(n, 1, 0.30)

  p_c <- exp(-1.35 + log(0.80) * treatment +
               0.006 * (age - 50) + 0.18 * comorbidity)
  p_d <- exp(-2.70 + log(0.60) * treatment +
               0.006 * (age - 50) + 0.18 * comorbidity)

  expect_true(all(p_d <= p_c))

  u <- runif(n)
  critical <- as.integer(u < p_d)
  composite <- as.integer(u < p_c)

  dat <- data.frame(
    treatment = treatment,
    age = age,
    comorbidity = comorbidity,
    composite = composite,
    critical = critical
  )

  fit_c <- glm(
    composite ~ treatment + age + comorbidity,
    family = poisson(link = "log"),
    data = dat
  )

  fit_d <- glm(
    critical ~ treatment + age + comorbidity,
    family = poisson(link = "log"),
    data = dat
  )

  x <- baco_models(
    fit_c,
    fit_d,
    exposure = "treatment",
    covariance_method = "joint_sandwich",
    effect_measure = "RR"
  )

  expect_true(x$nested_verified)
  expect_true(is.finite(x$primary$baco$estimate))
  expect_true(is.finite(x$primary$joint$covariance))
  expect_true(is.finite(x$primary$discrepancy$p_value))
})
