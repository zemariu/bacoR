#' BACO analysis from two correlated model coefficients
#'
#' Computes a BACO point estimate, Fieller confidence set, and the
#' composite-critical log-effect discrepancy from two correlated regression
#' coefficients estimated on the same effect scale.
#'
#' The sampling covariance between the coefficients is required. It cannot be
#' recovered from the two standard errors alone.
#'
#' @param beta_composite Estimated coefficient for the composite endpoint.
#'   The coefficient must be on a logarithmic effect scale (for example,
#'   log-RR, log-rate ratio, or log-HR).
#' @param beta_critical Estimated coefficient for the critical endpoint on the
#'   same effect scale as `beta_composite`.
#' @param se_composite Standard error of `beta_composite`.
#' @param se_critical Standard error of `beta_critical`.
#' @param covariance Sampling covariance between `beta_composite` and
#'   `beta_critical`.
#' @param conf_level Confidence level for Fieller and Delta intervals.
#' @param effect_measure Character label for the exponentiated effect measure,
#'   e.g. `"RR"`, `"IRR"`, or `"HR"`. This is used only for printing.
#' @param interpret Logical; if `TRUE`, print a plain-language interpretation
#'   of BACO, Delta, Fieller inference, and the joint discrepancy test.
#'
#' @details
#' The adjusted or model-based BACO is
#' `beta_composite / beta_critical`.
#'
#' The null hypothesis corresponding to BACO = 1, whenever the ratio is
#' defined, is tested as
#' `H0: beta_composite - beta_critical = 0`
#' using the joint sampling covariance:
#'
#' `Var(Delta) = Var(beta_composite) + Var(beta_critical) - 2*Cov`.
#'
#' This is not the conventional Wald test of the BACO ratio.
#'
#' @return An object of class `baco_coef_result`.
#' @export
baco_coef <- function(beta_composite,
                      beta_critical,
                      se_composite,
                      se_critical,
                      covariance,
                      conf_level = 0.95,
                      effect_measure = "effect",
                      interpret = FALSE) {

  args <- list(
    beta_composite = beta_composite,
    beta_critical = beta_critical,
    se_composite = se_composite,
    se_critical = se_critical,
    covariance = covariance,
    conf_level = conf_level
  )

  if (any(lengths(args) != 1L)) {
    stop("All numerical arguments must be scalar.", call. = FALSE)
  }

  vals <- unlist(args, use.names = FALSE)
  if (any(!is.finite(vals))) {
    stop("All numerical arguments must be finite.", call. = FALSE)
  }

  if (se_composite < 0 || se_critical < 0) {
    stop("Standard errors must be non-negative.", call. = FALSE)
  }

  if (conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must lie strictly between 0 and 1.", call. = FALSE)
  }

  if (length(effect_measure) != 1L || is.na(effect_measure) ||
      !nzchar(as.character(effect_measure))) {
    stop("`effect_measure` must be one non-empty character value.",
         call. = FALSE)
  }
  effect_measure <- as.character(effect_measure)

  var_composite <- se_composite^2
  var_critical <- se_critical^2

  max_cov <- sqrt(var_composite * var_critical)
  tol <- 100 * sqrt(.Machine$double.eps)

  if (abs(covariance) > max_cov + tol) {
    stop(
      paste0(
        "The supplied covariance is incompatible with the two standard errors. ",
        "|covariance| must not exceed se_composite * se_critical."
      ),
      call. = FALSE
    )
  }

  correlation <- if (se_composite > 0 && se_critical > 0) {
    covariance / (se_composite * se_critical)
  } else {
    NA_real_
  }

  effect_composite <- exp(beta_composite)
  effect_critical <- exp(beta_critical)

  baco_estimate <- if (abs(beta_critical) > sqrt(.Machine$double.eps)) {
    beta_composite / beta_critical
  } else {
    NA_real_
  }

  fieller <- baco_fieller(
    beta_composite = beta_composite,
    beta_critical = beta_critical,
    var_composite = var_composite,
    var_critical = var_critical,
    covariance = covariance,
    conf_level = conf_level
  )

  delta <- beta_composite - beta_critical
  var_delta <- var_composite + var_critical - 2 * covariance

  if (var_delta < 0 && abs(var_delta) < 1e-12) {
    var_delta <- 0
  }
  if (var_delta < 0) {
    stop("The implied variance of Delta is negative.", call. = FALSE)
  }

  se_delta <- sqrt(var_delta)
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)

  delta_ci <- c(
    lower = delta - zcrit * se_delta,
    upper = delta + zcrit * se_delta
  )

  delta_z <- if (se_delta > 0) delta / se_delta else NA_real_
  delta_p <- if (is.finite(delta_z)) {
    2 * stats::pnorm(-abs(delta_z))
  } else if (delta == 0) {
    1
  } else {
    NA_real_
  }

  # exp(Delta) = exponentiated composite effect / exponentiated critical effect.
  effect_ratio <- exp(delta)
  effect_ratio_ci <- exp(delta_ci)

  diagnostics <- list(
    baco_estimable = is.finite(baco_estimate),
    weak_denominator = fieller$weak_denominator,
    fieller_set_type = fieller$set_type,
    denominator_ci = fieller$denominator_ci,
    covariance_supplied = TRUE
  )

  out <- structure(
    list(
      call = match.call(),
      effect_measure = effect_measure,
      coefficients = list(
        beta_composite = beta_composite,
        beta_critical = beta_critical,
        se_composite = se_composite,
        se_critical = se_critical
      ),
      effects = list(
        composite = effect_composite,
        critical = effect_critical
      ),
      joint = list(
        covariance = covariance,
        correlation = correlation,
        var_composite = var_composite,
        var_critical = var_critical
      ),
      baco = list(
        estimate = baco_estimate,
        fieller = fieller
      ),
      discrepancy = list(
        delta = delta,
        se_delta = se_delta,
        ci = delta_ci,
        z = delta_z,
        p_value = delta_p,
        effect_ratio = effect_ratio,
        effect_ratio_ci = effect_ratio_ci
      ),
      diagnostics = diagnostics,
      interpret_enabled = isTRUE(interpret)
    ),
    class = "baco_coef_result"
  )

  out$interpretation <- .baco_build_model_interpretation(
    result = out,
    alpha = 1 - conf_level
  )

  out
}
