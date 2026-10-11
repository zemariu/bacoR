#' @export
print.baco_coef_result <- function(x, digits = 3,
                                   interpret = x$interpret_enabled, ...) {
  cat("BACO analysis from correlated model coefficients\n\n")

  measure <- x$effect_measure
  f <- x$baco$fieller

  cat("Model-based effects\n")
  cat(
    "  Composite ", measure, "                  : ",
    .baco_fmt_num(x$effects$composite, digits),
    "\n",
    sep = ""
  )
  cat(
    "  Critical ", measure, "                   : ",
    .baco_fmt_num(x$effects$critical, digits),
    "\n",
    sep = ""
  )

  cat("\nBACO ratio\n")

  if (identical(f$set_type, "bounded") && is.finite(x$baco$estimate)) {
    cat(
      "  Estimate [",
      round(100 * f$conf_level),
      "% Fieller CI]     : ",
      .baco_fmt_num(x$baco$estimate, digits),
      " ",
      .baco_format_fieller_set(f, digits),
      "\n",
      sep = ""
    )
  } else {
    cat(
      "  Estimate                     : ",
      if (is.finite(x$baco$estimate))
        .baco_fmt_num(x$baco$estimate, digits)
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
    .baco_fmt_num(x$discrepancy$delta, digits),
    " [",
    .baco_fmt_num(x$discrepancy$ci["lower"], digits),
    ", ",
    .baco_fmt_num(x$discrepancy$ci["upper"], digits),
    "]\n",
    sep = ""
  )
  cat(
    "  ", measure, " ratio (C / D)              : ",
    .baco_fmt_num(x$discrepancy$effect_ratio, digits),
    "\n",
    sep = ""
  )
  cat(
    "  Joint test p                  : ",
    if (is.finite(x$discrepancy$p_value))
      format.pval(x$discrepancy$p_value, digits = digits)
    else
      "not available",
    "\n",
    sep = ""
  )

  cat("\nJoint dependence\n")
  cat(
    "  Covariance                    : ",
    formatC(x$joint$covariance, digits = digits + 2, format = "f"),
    "\n",
    sep = ""
  )
  cat(
    "  Correlation                   : ",
    if (is.finite(x$joint$correlation))
      .baco_fmt_num(x$joint$correlation, digits)
    else
      "not available",
    "\n",
    sep = ""
  )

  cat("\nDiagnostics\n")
  cat("  Fieller set type              : ",
      x$diagnostics$fieller_set_type, "\n", sep = "")
  cat("  Weak denominator              : ",
      x$diagnostics$weak_denominator, "\n", sep = "")

  if (isTRUE(x$diagnostics$weak_denominator)) {
    cat(
      "  Warning                       : critical-effect coefficient is compatible with zero;\n",
      "                                  the BACO ratio is weakly determined.\n",
      sep = ""
    )
  }

  cat("\n")
  cat("Note: the p-value tests equality of the two model coefficients\n")
  cat("      using their joint covariance; it is not the conventional\n")
  cat("      Wald test of the BACO ratio.\n")


  if (isTRUE(interpret)) {
    .baco_print_interpretation(x, include_qrr = FALSE)
  }

  invisible(x)
}
