#' bacoR: Robust Inference for Bias Attributable to Composite Outcomes
#'
#' Tools for estimating and interpreting the bias attributable to composite
#' outcomes (BACO) index while distinguishing ratio estimation from
#' composite-critical effect discrepancy testing.
#'
#' The package provides three main interfaces:
#'
#' * `baco_counts()` for nested binary aggregate counts;
#' * `baco_effects()` and `baco_coef()` for externally estimated effects or
#'   regression coefficients with known joint covariance;
#' * `baco_models()` for two fitted log-link GLMs using the same participants.
#'
#' BACO is defined as the ratio of the logarithm of the association measure for
#' the composite endpoint to the logarithm of the association measure for the
#' critical component. Fieller confidence sets are used for the BACO ratio.
#' The no-discrepancy hypothesis is tested directly rather than through the
#' conventional Wald ratio test.
#'
#' Plain-language interpretation is optional and is disabled by default
#' (`interpret = FALSE`) in the public analytical interfaces.
#'
#' @keywords internal
#' @name bacoR-package
#' @aliases bacoR
"_PACKAGE"
