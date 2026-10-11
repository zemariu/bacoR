#' @export
print.baco_models_result <- function(x, digits = 3,
                                     interpret = x$interpret_enabled, ...) {
  cat("BACO analysis from fitted regression models\n\n")

  cat("Model information\n")
  cat("  Observations                   : ", x$n_observations, "\n", sep = "")
  cat("  Exposure coefficient          : ", x$coefficient_composite, "\n", sep = "")
  cat("  Covariance method             : ", x$covariance_method, "\n", sep = "")
  if (!is.na(x$nested_verified)) {
    cat("  Nested endpoints              : ", x$nested_verified, "\n", sep = "")
  }

  z <- x$primary
  f <- z$baco$fieller
  measure <- x$effect_measure

  cat("\nAdjusted/model-based effects\n")
  cat(
    "  Composite ", measure, "                  : ",
    .baco_fmt_num(z$effects$composite, digits),
    "\n",
    sep = ""
  )
  cat(
    "  Critical ", measure, "                   : ",
    .baco_fmt_num(z$effects$critical, digits),
    "\n",
    sep = ""
  )

  cat("\nBACO ratio\n")

  if (identical(f$set_type, "bounded") && is.finite(z$baco$estimate)) {
    cat(
      "  Estimate [",
      round(100 * f$conf_level),
      "% Fieller CI]     : ",
      .baco_fmt_num(z$baco$estimate, digits),
      " ",
      .baco_format_fieller_set(f, digits),
      "\n",
      sep = ""
    )
  } else {
    cat(
      "  Estimate                     : ",
      if (is.finite(z$baco$estimate))
        .baco_fmt_num(z$baco$estimate, digits)
      else
        "not estimable",
      "\n",
      sep = ""
    )
    cat(
      "  ",
      round(100 * f$conf_level),
      "% Fieller confidence set   : ",
      .baco_format_fieller_set(f, digits),
      "\n",
      sep = ""
    )
  }

  cat("\nComposite-critical discrepancy\n")
  cat(
    "  Delta [",
    round(100 * f$conf_level),
    "% CI]             : ",
    .baco_fmt_num(z$discrepancy$delta, digits),
    " [",
    .baco_fmt_num(z$discrepancy$ci["lower"], digits),
    ", ",
    .baco_fmt_num(z$discrepancy$ci["upper"], digits),
    "]\n",
    sep = ""
  )
  cat(
    "  ", measure, " ratio (C / D)              : ",
    .baco_fmt_num(z$discrepancy$effect_ratio, digits),
    "\n",
    sep = ""
  )
  cat(
    "  Joint test p                  : ",
    if (is.finite(z$discrepancy$p_value))
      format.pval(z$discrepancy$p_value, digits = digits)
    else
      "not available",
    "\n",
    sep = ""
  )

  cat("\nJoint dependence\n")
  cat(
    "  Covariance                    : ",
    formatC(z$joint$covariance, digits = digits + 2, format = "f"),
    "\n",
    sep = ""
  )
  cat(
    "  Correlation                   : ",
    if (is.finite(z$joint$correlation))
      .baco_fmt_num(z$joint$correlation, digits)
    else
      "not available",
    "\n",
    sep = ""
  )

  if (!is.null(x$bootstrap_info)) {
    cat("\nBootstrap validation\n")
    cat(
      "  Successful refits             : ",
      x$bootstrap_info$n_success,
      "/",
      x$bootstrap_info$R,
      " (",
      .baco_fmt_num(100 * x$bootstrap_info$success_fraction, 1),
      "%)\n",
      sep = ""
    )
    cat(
      "  Bootstrap covariance          : ",
      formatC(x$bootstrap_info$covariance, digits = digits + 2, format = "f"),
      "\n",
      sep = ""
    )
    cat(
      "  Bootstrap correlation         : ",
      .baco_fmt_num(x$bootstrap_info$correlation, digits),
      "\n",
      sep = ""
    )
  }

  cat("\nDiagnostics\n")
  cat("  Fieller set type              : ", f$set_type, "\n", sep = "")
  cat("  Weak denominator              : ", f$weak_denominator, "\n", sep = "")
  if (!is.na(x$nested_verified)) {
    cat("  Nested endpoints              : ", x$nested_verified, "\n", sep = "")
  }

  if (isTRUE(f$weak_denominator)) {
    cat(
      "  Warning                       : critical-effect coefficient is compatible with zero;\n",
      "                                  the BACO ratio is weakly determined.\n",
      sep = ""
    )
  }

  cat("\n")
  cat("Note: the reported p-value tests equality of the two exposure\n")
  cat("      coefficients using their joint covariance. It is not the\n")
  cat("      conventional Wald test of the BACO ratio.\n")


  if (isTRUE(interpret)) {
    .baco_print_interpretation(z, include_qrr = FALSE)
  }

  invisible(x)
}
