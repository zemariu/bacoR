.baco_fmt_num <- function(x, digits = 3) {
  if (!is.finite(x)) {
    if (is.infinite(x) && x < 0) return("-Inf")
    if (is.infinite(x) && x > 0) return("Inf")
    return("NA")
  }
  formatC(x, digits = digits, format = "f")
}

.baco_format_fieller_set <- function(f, digits = 3) {
  if (is.null(f) || is.null(f$set_type)) {
    return("not available")
  }

  if (identical(f$set_type, "empty")) {
    return("empty")
  }

  ints <- f$intervals
  if (is.null(ints) || nrow(ints) == 0L) {
    return(f$set_type)
  }

  pieces <- character(nrow(ints))

  for (i in seq_len(nrow(ints))) {
    lo <- ints$lower[i]
    hi <- ints$upper[i]

    left_bracket <- if (is.infinite(lo)) "(" else "["
    right_bracket <- if (is.infinite(hi)) ")" else "]"

    pieces[i] <- paste0(
      left_bracket,
      .baco_fmt_num(lo, digits),
      ", ",
      .baco_fmt_num(hi, digits),
      right_bracket
    )
  }

  paste(pieces, collapse = " U ")
}


.baco_print_interpretation <- function(x, include_qrr = FALSE) {
  if (is.null(x$interpretation)) return(invisible(NULL))

  cat("\nInterpretation\n")
  cat("  BACO                          : ",
      x$interpretation$baco$label, "\n", sep = "")
  cat("                                  ",
      x$interpretation$baco$text, "\n", sep = "")

  cat("  Discrepancy                   : ",
      x$interpretation$discrepancy$text, "\n", sep = "")

  if (isTRUE(include_qrr) && !is.null(x$interpretation$qrr)) {
    cat("  QRR                           : ",
        x$interpretation$qrr$text, "\n", sep = "")
  }

  cat("  Testing                       : ",
      x$interpretation$testing$text, "\n", sep = "")

  if (!is.null(x$interpretation$fieller)) {
    cat("  Fieller                       : ",
        x$interpretation$fieller$text, "\n", sep = "")
  }

  invisible(NULL)
}
