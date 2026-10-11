# Internal helper for a probability-ordered two-sided Fisher mid-P value.
.baco_midp <- function(d_t, c_t, d_c, c_c) {
  if (c_t <= 0 || c_c <= 0) return(NA_real_)

  d_total <- d_t + d_c
  if (d_total == 0) return(1)

  lo <- max(0, d_total - c_c)
  hi <- min(c_t, d_total)
  support <- lo:hi

  probs <- stats::dhyper(support, m = c_t, n = c_c, k = d_total)
  p_obs <- stats::dhyper(d_t, m = c_t, n = c_c, k = d_total)

  p_exact <- sum(probs[probs <= p_obs * (1 + 1e-12)])
  min(1, max(0, p_exact - 0.5 * p_obs))
}

.baco_scalar_count <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x < 0 || abs(x - round(x)) > 1e-8) {
    stop(sprintf("`%s` must be one non-negative integer-like value.", name),
         call. = FALSE)
  }
  as.numeric(x)
}

#' BACO analysis from nested binary aggregate counts
#'
#' Computes the conventional BACO point estimate when defined, Fieller
#' confidence sets for the BACO ratio, the log-effect discrepancy Delta,
#' QRR, and an exact conditional Fisher test of no composite-critical effect
#' discrepancy.
#'
#' The critical endpoint must be nested within the composite endpoint.
#' No continuity correction is applied.
#'
#' @param n_treatment Number randomized/observed in the treatment group.
#' @param n_control Number randomized/observed in the control group.
#' @param composite_treatment Composite endpoint events in treatment.
#' @param composite_control Composite endpoint events in control.
#' @param critical_treatment Critical endpoint events in treatment.
#' @param critical_control Critical endpoint events in control.
#' @param conf_level Confidence level for Fieller and Delta intervals.
#' @param midp Logical; also calculate the probability-ordered Fisher mid-P.
#' @param interpret Logical; if `TRUE`, print a plain-language interpretation
#'   of BACO, Delta/QRR, Fieller inference, and the exact discrepancy test.
#'
#' @return An object of class `baco_result`.
#' @export
baco_counts <- function(n_treatment,
                        n_control,
                        composite_treatment,
                        composite_control,
                        critical_treatment,
                        critical_control,
                        conf_level = 0.95,
                        midp = FALSE,
                        interpret = FALSE) {

  n_treatment <- .baco_scalar_count(n_treatment, "n_treatment")
  n_control <- .baco_scalar_count(n_control, "n_control")
  composite_treatment <- .baco_scalar_count(composite_treatment, "composite_treatment")
  composite_control <- .baco_scalar_count(composite_control, "composite_control")
  critical_treatment <- .baco_scalar_count(critical_treatment, "critical_treatment")
  critical_control <- .baco_scalar_count(critical_control, "critical_control")

  if (n_treatment <= 0 || n_control <= 0) {
    stop("Group sample sizes must be greater than zero.", call. = FALSE)
  }
  if (composite_treatment > n_treatment || composite_control > n_control) {
    stop("Composite-event counts cannot exceed group sample sizes.", call. = FALSE)
  }
  if (critical_treatment > composite_treatment || critical_control > composite_control) {
    stop("Critical events must be nested within composite events (D <= C).",
         call. = FALSE)
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must lie strictly between 0 and 1.", call. = FALSE)
  }

  risk_c_t <- composite_treatment / n_treatment
  risk_c_c <- composite_control / n_control
  risk_d_t <- critical_treatment / n_treatment
  risk_d_c <- critical_control / n_control

  rr_composite <- if (risk_c_c == 0) {
    if (risk_c_t == 0) NA_real_ else Inf
  } else {
    risk_c_t / risk_c_c
  }

  rr_critical <- if (risk_d_c == 0) {
    if (risk_d_t == 0) NA_real_ else Inf
  } else {
    risk_d_t / risk_d_c
  }

  beta_composite <- if (is.finite(rr_composite) && rr_composite > 0) {
    log(rr_composite)
  } else {
    NA_real_
  }

  beta_critical <- if (is.finite(rr_critical) && rr_critical > 0) {
    log(rr_critical)
  } else {
    NA_real_
  }

  var_composite <- if (composite_treatment > 0 && composite_control > 0) {
    1 / composite_treatment + 1 / composite_control -
      1 / n_treatment - 1 / n_control
  } else {
    NA_real_
  }

  var_critical <- if (critical_treatment > 0 && critical_control > 0) {
    1 / critical_treatment + 1 / critical_control -
      1 / n_treatment - 1 / n_control
  } else {
    NA_real_
  }

  # Nested binary aggregate-outcome covariance.
  covariance <- if (is.finite(var_composite)) var_composite else NA_real_

  baco_estimate <- if (is.finite(beta_composite) &&
                       is.finite(beta_critical) &&
                       abs(beta_critical) > sqrt(.Machine$double.eps)) {
    beta_composite / beta_critical
  } else {
    NA_real_
  }

  fieller <- NULL
  if (all(is.finite(c(beta_composite, beta_critical,
                      var_composite, var_critical, covariance)))) {
    fieller <- baco_fieller(
      beta_composite = beta_composite,
      beta_critical = beta_critical,
      var_composite = var_composite,
      var_critical = var_critical,
      covariance = covariance,
      conf_level = conf_level
    )
  }

  delta <- if (is.finite(beta_composite) && is.finite(beta_critical)) {
    beta_composite - beta_critical
  } else {
    NA_real_
  }

  var_delta <- if (all(is.finite(c(var_composite, var_critical, covariance)))) {
    var_composite + var_critical - 2 * covariance
  } else {
    NA_real_
  }

  if (is.finite(var_delta) && var_delta < 0 && abs(var_delta) < 1e-12) {
    var_delta <- 0
  }

  se_delta <- if (is.finite(var_delta) && var_delta > 0) sqrt(var_delta) else NA_real_
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  delta_ci <- if (is.finite(delta) && is.finite(se_delta)) {
    c(lower = delta - z * se_delta, upper = delta + z * se_delta)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }

  delta_p_wald <- if (is.finite(delta) && is.finite(se_delta) && se_delta > 0) {
    2 * stats::pnorm(-abs(delta / se_delta))
  } else {
    NA_real_
  }

  qrr <- if (composite_treatment > 0 && composite_control > 0) {
    q_t <- critical_treatment / composite_treatment
    q_c <- critical_control / composite_control
    if (q_c == 0) {
      if (q_t == 0) NA_real_ else Inf
    } else {
      q_t / q_c
    }
  } else {
    NA_real_
  }

  fisher_p <- NA_real_
  midp_p <- NA_real_
  fisher_note <- NULL

  if (composite_treatment > 0 && composite_control > 0) {
    if (critical_treatment == 0 && critical_control == 0) {
      fisher_p <- 1
      midp_p <- if (isTRUE(midp)) 1 else NA_real_
      fisher_note <- paste(
        "No critical events occurred in either group; p = 1 reflects",
        "absence of critical-event information, not evidence of equivalence."
      )
    } else {
      tab <- matrix(
        c(
          critical_treatment,
          composite_treatment - critical_treatment,
          critical_control,
          composite_control - critical_control
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
          group = c("treatment", "control"),
          event = c("critical", "noncritical_composite")
        )
      )
      fisher_p <- stats::fisher.test(tab, alternative = "two.sided")$p.value
      if (isTRUE(midp)) {
        midp_p <- .baco_midp(
          critical_treatment,
          composite_treatment,
          critical_control,
          composite_control
        )
      }
    }
  } else {
    fisher_note <- paste(
      "Exact discrepancy testing is not informative because at least one",
      "group had no composite events."
    )
  }

  diagnostics <- list(
    total_critical_events = critical_treatment + critical_control,
    zero_critical_any_arm = critical_treatment == 0 || critical_control == 0,
    zero_critical_both_arms = critical_treatment == 0 && critical_control == 0,
    zero_composite_any_arm = composite_treatment == 0 || composite_control == 0,
    baco_estimable = is.finite(baco_estimate),
    fieller_available = !is.null(fieller),
    fieller_set_type = if (!is.null(fieller)) fieller$set_type else NA_character_,
    weak_denominator = if (!is.null(fieller)) fieller$weak_denominator else NA
  )

  interpretation <- list(
    baco = .baco_original_interpretation(
      theta = baco_estimate,
      weak_denominator = diagnostics$weak_denominator
    ),
    discrepancy = if (is.finite(delta) && is.finite(rr_composite) &&
                      is.finite(rr_critical)) {
      .baco_effect_discrepancy_interpretation(
        effect_composite = rr_composite,
        effect_critical = rr_critical,
        delta = delta,
        effect_ratio = exp(delta),
        effect_measure = "RR"
      )
    } else {
      list(
        label = "Not estimable",
        text = "Delta is not estimable from the observed risk ratios."
      )
    },
    qrr = .baco_qrr_interpretation(qrr),
    testing = .baco_test_interpretation(
      p_value = fisher_p,
      alpha = 1 - conf_level,
      method = "Fisher exact test"
    ),
    fieller = .baco_fieller_interpretation(fieller)
  )

  structure(
    list(
      call = match.call(),
      input = list(
        n_treatment = n_treatment,
        n_control = n_control,
        composite_treatment = composite_treatment,
        composite_control = composite_control,
        critical_treatment = critical_treatment,
        critical_control = critical_control
      ),
      effects = list(
        rr_composite = rr_composite,
        rr_critical = rr_critical,
        beta_composite = beta_composite,
        beta_critical = beta_critical
      ),
      baco = list(
        estimate = baco_estimate,
        fieller = fieller
      ),
      discrepancy = list(
        delta = delta,
        se_delta = se_delta,
        delta_ci = delta_ci,
        delta_p_wald = delta_p_wald,
        qrr = qrr
      ),
      exact_test = list(
        fisher_p = fisher_p,
        midp_p = midp_p,
        note = fisher_note
      ),
      covariance = list(
        var_composite = var_composite,
        var_critical = var_critical,
        covariance = covariance
      ),
      diagnostics = diagnostics,
      interpretation = interpretation,
      interpret_enabled = isTRUE(interpret)
    ),
    class = "baco_result"
  )
}
