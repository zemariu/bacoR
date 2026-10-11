#' @export
print.baco_result <- function(x, digits = 3,
                              interpret = x$interpret_enabled, ...) {
  cat("BACO analysis from nested aggregate counts\n\n")

  cat("Effects\n")
  cat(
    "  Composite RR                 : ",
    .baco_fmt_num(x$effects$rr_composite, digits),
    "\n",
    sep = ""
  )
  cat(
    "  Critical RR                  : ",
    .baco_fmt_num(x$effects$rr_critical, digits),
    "\n",
    sep = ""
  )

  cat("\nBACO ratio\n")

  if (!is.null(x$baco$fieller) &&
      identical(x$baco$fieller$set_type, "bounded") &&
      is.finite(x$baco$estimate)) {

    f <- x$baco$fieller
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

    if (!is.null(x$baco$fieller)) {
      f <- x$baco$fieller
      cat(
        "  ",
        round(100 * f$conf_level),
        "% Fieller confidence set   : ",
        .baco_format_fieller_set(f, digits),
        "\n",
        sep = ""
      )
    } else {
      cat(
        "  95% Fieller confidence set   : not available from observed counts\n"
      )
    }
  }

  cat("\nComposite-critical discrepancy\n")
  cat(
    "  Delta                         : ",
    if (is.finite(x$discrepancy$delta))
      .baco_fmt_num(x$discrepancy$delta, digits)
    else
      "not estimable",
    "\n",
    sep = ""
  )
  cat(
    "  QRR                           : ",
    if (!is.na(x$discrepancy$qrr))
      .baco_fmt_num(x$discrepancy$qrr, digits)
    else
      "not estimable",
    "\n",
    sep = ""
  )

  cat("\nExact null-discrepancy test\n")
  cat(
    "  Fisher exact p                : ",
    if (is.finite(x$exact_test$fisher_p))
      format.pval(x$exact_test$fisher_p, digits = digits)
    else
      "not available",
    "\n",
    sep = ""
  )

  if (is.finite(x$exact_test$midp_p)) {
    cat(
      "  Mid-P (sensitivity)           : ",
      format.pval(x$exact_test$midp_p, digits = digits),
      "\n",
      sep = ""
    )
  }

  if (!is.null(x$exact_test$note)) {
    cat("  Note                          : ", x$exact_test$note, "\n", sep = "")
  }

  cat("\nDiagnostics\n")
  if (!is.null(x$baco$fieller)) {
    cat("  Fieller set type              : ",
        x$baco$fieller$set_type, "\n", sep = "")
  } else {
    cat("  Fieller set type              : not available\n")
  }
  cat("  Weak denominator              : ",
      x$diagnostics$weak_denominator, "\n", sep = "")
  cat("  Total critical events         : ",
      x$diagnostics$total_critical_events, "\n", sep = "")
  cat("  Zero critical cell            : ",
      x$diagnostics$zero_critical_any_arm, "\n", sep = "")


  if (isTRUE(interpret)) {
    .baco_print_interpretation(x, include_qrr = TRUE)
  }

  invisible(x)
}
