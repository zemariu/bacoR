#' Fieller confidence set for the BACO ratio
#'
#' Computes a Fieller confidence set for the ratio
#' `beta_composite / beta_critical`, explicitly incorporating the covariance
#' between the two effect estimates.
#'
#' @param beta_composite Estimated log-effect for the composite endpoint.
#' @param beta_critical Estimated log-effect for the critical endpoint.
#' @param var_composite Sampling variance of `beta_composite`.
#' @param var_critical Sampling variance of `beta_critical`.
#' @param covariance Sampling covariance between the two effect estimates.
#' @param conf_level Confidence level. Default is 0.95.
#' @param tol Numerical tolerance.
#'
#' @return A list containing the BACO point estimate, Fieller set
#'   classification, interval components, denominator confidence interval,
#'   and weak-denominator diagnostic.
#' @export
baco_fieller <- function(beta_composite,
                         beta_critical,
                         var_composite,
                         var_critical,
                         covariance,
                         conf_level = 0.95,
                         tol = sqrt(.Machine$double.eps)) {

  vals <- c(beta_composite, beta_critical, var_composite,
            var_critical, covariance, conf_level)

  if (any(lengths(list(beta_composite, beta_critical, var_composite,
                       var_critical, covariance, conf_level)) != 1L)) {
    stop("All arguments must be scalar.", call. = FALSE)
  }

  if (any(!is.finite(vals))) {
    stop("All arguments must be finite.", call. = FALSE)
  }

  if (var_composite < 0 || var_critical < 0) {
    stop("Variances must be non-negative.", call. = FALSE)
  }

  if (conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must lie strictly between 0 and 1.", call. = FALSE)
  }

  # Numerical check of positive semidefiniteness of the 2x2 covariance matrix.
  max_cov <- sqrt(var_composite * var_critical)
  if (abs(covariance) > max_cov + 100 * tol) {
    stop("The supplied covariance matrix is not positive semidefinite.",
         call. = FALSE)
  }

  alpha <- 1 - conf_level
  z <- stats::qnorm(1 - alpha / 2)
  z2 <- z^2

  A <- beta_critical^2 - z2 * var_critical
  B <- -2 * beta_composite * beta_critical + 2 * z2 * covariance
  C <- beta_composite^2 - z2 * var_composite

  discriminant <- B^2 - 4 * A * C
  disc_tol <- tol * max(1, abs(B)^2, abs(4 * A * C))

  point_estimate <- if (abs(beta_critical) > tol) {
    beta_composite / beta_critical
  } else {
    NA_real_
  }

  se_critical <- sqrt(var_critical)
  denominator_ci <- c(
    lower = beta_critical - z * se_critical,
    upper = beta_critical + z * se_critical
  )

  # A <= 0 is the usual Fieller signal that zero is not excluded for the
  # denominator effect at the requested confidence level.
  weak_denominator <- A <= tol

  intervals <- data.frame(lower = numeric(0), upper = numeric(0))
  set_type <- NA_character_

  if (abs(A) <= tol) {
    # Degenerate quadratic -> linear inequality B*theta + C <= 0.
    if (abs(B) <= tol) {
      if (C <= tol) {
        set_type <- "all_real"
        intervals <- data.frame(lower = -Inf, upper = Inf)
      } else {
        set_type <- "empty"
      }
    } else {
      root <- -C / B
      set_type <- "half_line"
      if (B > 0) {
        intervals <- data.frame(lower = -Inf, upper = root)
      } else {
        intervals <- data.frame(lower = root, upper = Inf)
      }
    }
  } else if (A > 0) {
    if (discriminant < -disc_tol) {
      set_type <- "empty"
    } else {
      discriminant <- max(discriminant, 0)
      root1 <- (-B - sqrt(discriminant)) / (2 * A)
      root2 <- (-B + sqrt(discriminant)) / (2 * A)
      roots <- sort(c(root1, root2))
      if (abs(diff(roots)) <= tol) {
        set_type <- "point"
      } else {
        set_type <- "bounded"
      }
      intervals <- data.frame(lower = roots[1], upper = roots[2])
    }
  } else {
    # A < 0: the quadratic opens downward.
    if (discriminant <= disc_tol) {
      set_type <- "all_real"
      intervals <- data.frame(lower = -Inf, upper = Inf)
    } else {
      root1 <- (-B - sqrt(discriminant)) / (2 * A)
      root2 <- (-B + sqrt(discriminant)) / (2 * A)
      roots <- sort(c(root1, root2))
      set_type <- "disjoint_unbounded"
      intervals <- data.frame(
        lower = c(-Inf, roots[2]),
        upper = c(roots[1], Inf)
      )
    }
  }

  structure(
    list(
      estimate = point_estimate,
      conf_level = conf_level,
      set_type = set_type,
      intervals = intervals,
      weak_denominator = weak_denominator,
      denominator_ci = denominator_ci,
      beta_composite = beta_composite,
      beta_critical = beta_critical,
      var_composite = var_composite,
      var_critical = var_critical,
      covariance = covariance,
      quadratic = c(A = A, B = B, C = C),
      discriminant = discriminant
    ),
    class = "baco_fieller"
  )
}
