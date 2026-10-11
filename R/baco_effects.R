#' BACO analysis from exponentiated adjusted effects
#'
#' Convenience interface for users who have adjusted relative effects such as
#' risk ratios (RR), incidence-rate ratios (IRR), or hazard ratios (HR), plus
#' the standard errors of their underlying log coefficients.
#'
#' @param effect_composite Exponentiated adjusted effect for the composite
#'   endpoint, e.g. an adjusted RR.
#' @param effect_critical Exponentiated adjusted effect for the critical
#'   endpoint, on the same scale.
#' @param se_log_composite Standard error of log(`effect_composite`), i.e. the
#'   regression-coefficient standard error.
#' @param se_log_critical Standard error of log(`effect_critical`).
#' @param covariance_log Optional sampling covariance between the two log
#'   coefficients.
#' @param correlation_log Optional sampling correlation between the two log
#'   coefficients. If supplied instead of `covariance_log`, the covariance is
#'   calculated internally.
#' @param conf_level Confidence level.
#' @param effect_measure Label such as `"RR"`, `"IRR"`, or `"HR"`.
#' @param interpret Logical; if `TRUE`, print plain-language interpretation.
#'
#' @details
#' The two standard errors alone do not identify the covariance. Supply either
#' `covariance_log` or `correlation_log`.
#'
#' @return An object of class `baco_coef_result`.
#' @export
baco_effects <- function(effect_composite,
                         effect_critical,
                         se_log_composite,
                         se_log_critical,
                         covariance_log = NULL,
                         correlation_log = NULL,
                         conf_level = 0.95,
                         effect_measure = "RR",
                         interpret = FALSE) {

  if (length(effect_composite) != 1L ||
      length(effect_critical) != 1L ||
      !is.finite(effect_composite) ||
      !is.finite(effect_critical) ||
      effect_composite <= 0 ||
      effect_critical <= 0) {
    stop("Both adjusted effects must be finite positive scalar values.",
         call. = FALSE)
  }

  if (length(se_log_composite) != 1L ||
      length(se_log_critical) != 1L ||
      !is.finite(se_log_composite) ||
      !is.finite(se_log_critical) ||
      se_log_composite < 0 ||
      se_log_critical < 0) {
    stop("Log-scale standard errors must be finite non-negative scalars.",
         call. = FALSE)
  }

  if (is.null(covariance_log) && is.null(correlation_log)) {
    stop(
      paste0(
        "Supply either `covariance_log` or `correlation_log`. ",
        "The covariance cannot be recovered from the two standard errors alone."
      ),
      call. = FALSE
    )
  }

  if (!is.null(correlation_log)) {
    if (length(correlation_log) != 1L ||
        !is.finite(correlation_log) ||
        correlation_log < -1 ||
        correlation_log > 1) {
      stop("`correlation_log` must lie between -1 and 1.", call. = FALSE)
    }

    cov_from_cor <- correlation_log * se_log_composite * se_log_critical

    if (is.null(covariance_log)) {
      covariance_log <- cov_from_cor
    } else {
      if (length(covariance_log) != 1L || !is.finite(covariance_log)) {
        stop("`covariance_log` must be one finite value.", call. = FALSE)
      }

      tol <- 1e-8 * max(1, abs(covariance_log), abs(cov_from_cor))
      if (abs(covariance_log - cov_from_cor) > tol) {
        stop(
          "`covariance_log` and `correlation_log` imply inconsistent dependence.",
          call. = FALSE
        )
      }
    }
  }

  ans <- baco_coef(
    beta_composite = log(effect_composite),
    beta_critical = log(effect_critical),
    se_composite = se_log_composite,
    se_critical = se_log_critical,
    covariance = covariance_log,
    conf_level = conf_level,
    effect_measure = effect_measure,
    interpret = interpret
  )

  ans$input_scale <- "exponentiated_effects"
  ans$input_effects <- c(
    composite = effect_composite,
    critical = effect_critical
  )
  ans$call <- match.call()
  ans
}
