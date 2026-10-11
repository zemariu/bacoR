#' Simulated nested-endpoint trial for BACO examples
#'
#' A deterministic simulated individual-level randomized trial intended for
#' examples, vignettes, unit tests, and demonstrations of `baco_models()`.
#'
#' The data-generating mechanism targets approximately:
#' \itemize{
#'   \item composite-endpoint RR = 0.80;
#'   \item critical-endpoint RR = 0.60;
#'   \item population BACO = \code{log(0.80) / log(0.60)};
#'   \item population Delta = \code{log(0.80) - log(0.60)}.
#' }
#'
#' The critical endpoint is nested within the composite endpoint by design,
#' such that \code{critical <= composite} for every participant.
#'
#' @format A data frame with 2,500 rows and 7 variables:
#' \describe{
#'   \item{id}{Participant identifier.}
#'   \item{treatment}{Randomized treatment indicator (0 = control, 1 = intervention).}
#'   \item{age}{Age in years.}
#'   \item{sex}{Binary sex indicator used as an adjustment covariate.}
#'   \item{comorbidity}{Binary comorbidity indicator.}
#'   \item{composite}{Composite endpoint indicator.}
#'   \item{critical}{Critical endpoint indicator nested within \code{composite}.}
#' }
#'
#' @details
#' The dataset is simulated and contains no real or confidential participant
#' information. It is included solely for documentation, examples, teaching,
#' vignettes, and software testing.
#'
#' The object also stores the simulation targets in the
#' \code{data_generating_parameters} attribute.
#'
#' @source Simulated specifically for the bacoR package.
#'
"baco_trial"
