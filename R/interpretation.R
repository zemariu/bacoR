# Internal interpretation helpers ---------------------------------------------

.baco_original_interpretation <- function(theta,
                                          weak_denominator = FALSE,
                                          tol = 1e-10) {

  if (!is.finite(theta)) {
    return(list(
      classification = "not_estimable",
      label = "Not estimable",
      text = paste(
        "The BACO ratio is not estimable, so the original BACO",
        "classification cannot be applied."
      )
    ))
  }

  if (abs(theta - 1) <= tol) {
    classification <- "no_bias"
    label <- "No bias attributable to the composite outcome"
    text <- paste(
      "According to the original Diaz-Quijano interpretation,",
      "BACO = 1 indicates no bias attributable to use of the composite outcome."
    )
  } else if (theta > 1) {
    classification <- "overestimation"
    label <- "Overestimation"
    text <- paste(
      "According to the original Diaz-Quijano interpretation,",
      "BACO > 1 indicates that the association with the composite outcome",
      "is stronger than that with the critical component."
    )
  } else if (theta > 0 && theta < 1) {
    classification <- "underestimation"
    label <- "Underestimation (toward the null)"
    text <- paste(
      "According to the original Diaz-Quijano interpretation,",
      "0 < BACO < 1 indicates that the composite-outcome association",
      "is weaker and biased toward the null relative to the critical component."
    )
  } else if (abs(theta) <= tol) {
    classification <- "underestimation_to_null"
    label <- "Composite association at the null"
    text <- paste(
      "BACO = 0 is a boundary case in which the composite association",
      "is at the null while the critical-component association is non-null;",
      "this is consistent with complete attenuation toward the null."
    )
  } else {
    classification <- "inversion"
    label <- "Inversion"
    text <- paste(
      "According to the original Diaz-Quijano interpretation,",
      "BACO < 0 indicates an inversion of the association relative to",
      "the critical component."
    )
  }

  if (isTRUE(weak_denominator)) {
    text <- paste(
      text,
      "However, the critical-effect denominator is compatible with zero;",
      "the BACO ratio is weakly determined and this point-estimate",
      "classification should be interpreted cautiously."
    )
  }

  list(
    classification = classification,
    label = label,
    text = text
  )
}


.baco_effect_discrepancy_interpretation <- function(effect_composite,
                                                    effect_critical,
                                                    delta,
                                                    effect_ratio,
                                                    effect_measure = "effect",
                                                    tol = 1e-10) {

  if (!all(is.finite(c(effect_composite, effect_critical, delta, effect_ratio)))) {
    return(list(
      label = "Not estimable",
      text = "The effect discrepancy is not estimable from the available data."
    ))
  }

  if (abs(delta) <= tol) {
    return(list(
      label = "No log-effect discrepancy",
      text = paste0(
        "Delta = 0 and the composite and critical ", effect_measure,
        " estimates are equal."
      )
    ))
  }

  pct <- abs(effect_ratio - 1) * 100

  if (effect_ratio > 1) {
    relative_text <- paste0(
      "The composite ", effect_measure, " is ",
      formatC(pct, digits = 1, format = "f"),
      "% higher than the critical ", effect_measure, "."
    )
  } else {
    relative_text <- paste0(
      "The composite ", effect_measure, " is ",
      formatC(pct, digits = 1, format = "f"),
      "% lower than the critical ", effect_measure, "."
    )
  }

  # Contextualize relative to the null without assuming causality.
  if (effect_composite < 1 && effect_critical < 1) {
    if (effect_composite > effect_critical) {
      null_text <- paste(
        "Because both estimates are below 1, the composite association",
        "is closer to the null (weaker in magnitude)."
      )
    } else {
      null_text <- paste(
        "Because both estimates are below 1, the composite association",
        "is farther from the null (stronger in magnitude)."
      )
    }
  } else if (effect_composite > 1 && effect_critical > 1) {
    if (effect_composite > effect_critical) {
      null_text <- paste(
        "Because both estimates are above 1, the composite association",
        "is farther from the null (stronger in magnitude)."
      )
    } else {
      null_text <- paste(
        "Because both estimates are above 1, the composite association",
        "is closer to the null (weaker in magnitude)."
      )
    }
  } else if ((effect_composite < 1 && effect_critical > 1) ||
             (effect_composite > 1 && effect_critical < 1)) {
    null_text <- paste(
      "The two association measures are on opposite sides of the null,",
      "indicating an inversion of direction."
    )
  } else {
    null_text <- ""
  }

  list(
    label = "Composite-critical effect discrepancy",
    text = trimws(paste(relative_text, null_text))
  )
}


.baco_qrr_interpretation <- function(qrr) {
  if (is.na(qrr)) {
    return(list(
      label = "QRR not estimable",
      text = paste(
        "The conditional critical-event risk ratio (QRR) is not estimable",
        "from the observed composite-event counts."
      )
    ))
  }

  if (is.infinite(qrr)) {
    return(list(
      label = "QRR infinite",
      text = paste(
        "Among participants with a composite event, the critical component",
        "was absent in the control group but present in the treatment group."
      )
    ))
  }

  if (abs(qrr - 1) <= 1e-10) {
    return(list(
      label = "Equal critical-event composition",
      text = paste(
        "Among participants with a composite event, the proportion with the",
        "critical component is the same in the two groups (QRR = 1)."
      )
    ))
  }

  pct <- abs(qrr - 1) * 100

  if (qrr < 1) {
    text <- paste0(
      "Among participants with a composite event, the proportion with the ",
      "critical component is ",
      formatC(pct, digits = 1, format = "f"),
      "% lower in the treatment group than in the control group."
    )
  } else {
    text <- paste0(
      "Among participants with a composite event, the proportion with the ",
      "critical component is ",
      formatC(pct, digits = 1, format = "f"),
      "% higher in the treatment group than in the control group."
    )
  }

  list(
    label = "Conditional critical-event composition",
    text = text
  )
}


.baco_test_interpretation <- function(p_value,
                                      alpha = 0.05,
                                      method = "test") {

  if (!is.finite(p_value)) {
    return(list(
      conclusion = "not_available",
      text = "The null-discrepancy test is not available for these data."
    ))
  }

  if (p_value < alpha) {
    return(list(
      conclusion = "evidence_of_discrepancy",
      text = paste0(
        "At alpha = ", formatC(alpha, digits = 3, format = "fg"),
        ", the ", method,
        " provides evidence that the composite and critical effects differ. ",
        "When BACO is defined, this tests the null hypothesis corresponding ",
        "to BACO = 1."
      )
    ))
  }

  list(
    conclusion = "no_evidence_of_discrepancy",
    text = paste0(
      "At alpha = ", formatC(alpha, digits = 3, format = "fg"),
      ", the ", method,
      " does not provide evidence that the composite and critical effects differ. ",
      "This does not establish equivalence or absence of bias."
    )
  )
}


.baco_fieller_interpretation <- function(fieller) {
  if (is.null(fieller)) {
    return(list(
      text = "A Fieller confidence set is not available."
    ))
  }

  if (identical(fieller$set_type, "bounded")) {
    return(list(
      text = paste(
        "The Fieller confidence set is bounded: the BACO ratio has a finite",
        "confidence interval and the critical-effect denominator excludes zero",
        "at the selected confidence level."
      )
    ))
  }

  if (isTRUE(fieller$weak_denominator)) {
    return(list(
      text = paste(
        "The critical-effect coefficient is compatible with zero.",
        "The BACO ratio is therefore weakly determined; an unbounded or",
        "nonstandard Fieller confidence set is an inferential result, not",
        "a computational failure."
      )
    ))
  }

  list(
    text = paste(
      "The Fieller confidence set is nonstandard and should be reported",
      "in its full form rather than forced into a finite interval."
    )
  )
}


.baco_build_model_interpretation <- function(result,
                                             alpha = 0.05) {

  list(
    baco = .baco_original_interpretation(
      theta = result$baco$estimate,
      weak_denominator = result$baco$fieller$weak_denominator
    ),
    discrepancy = .baco_effect_discrepancy_interpretation(
      effect_composite = result$effects$composite,
      effect_critical = result$effects$critical,
      delta = result$discrepancy$delta,
      effect_ratio = result$discrepancy$effect_ratio,
      effect_measure = result$effect_measure
    ),
    testing = .baco_test_interpretation(
      p_value = result$discrepancy$p_value,
      alpha = alpha,
      method = "joint coefficient test"
    ),
    fieller = .baco_fieller_interpretation(
      result$baco$fieller
    )
  )
}
