# Internal: resolve an exact coefficient or an unambiguous factor coefficient.
.baco_resolve_coef <- function(model, exposure) {
  nm <- names(stats::coef(model))

  if (exposure %in% nm) {
    return(exposure)
  }

  candidates <- nm[startsWith(nm, exposure)]
  if (length(candidates) == 1L) {
    return(candidates)
  }

  stop(
    paste0(
      "Could not identify a unique coefficient for `", exposure, "`. ",
      "Available coefficients are: ", paste(nm, collapse = ", "), ". ",
      "Supply the exact coefficient name."
    ),
    call. = FALSE
  )
}

# Internal: individual score contributions and HC0 sandwich ingredients for GLM.
.baco_glm_sandwich_parts <- function(model) {
  if (!inherits(model, "glm")) {
    stop("The current `baco_models()` implementation supports `glm` objects only.",
         call. = FALSE)
  }

  mf <- stats::model.frame(model)
  y <- stats::model.response(mf)

  if (!is.numeric(y) || is.matrix(y)) {
    stop(
      "The current model interface requires a numeric individual-level response.",
      call. = FALSE
    )
  }

  X <- stats::model.matrix(model)
  mu <- stats::fitted(model)
  eta <- stats::predict(model, type = "link")

  fam <- model$family
  dmu <- fam$mu.eta(eta)
  vv <- fam$variance(mu)

  pw <- model$prior.weights
  if (is.null(pw)) pw <- rep(1, length(mu))

  if (any(!is.finite(c(mu, eta, dmu, vv, pw))) || any(vv <= 0)) {
    stop("Non-finite GLM quantities prevented joint covariance estimation.",
         call. = FALSE)
  }

  # Quasi-score contribution. Dispersion cancels from the sandwich expression,
  # so it is intentionally omitted here.
  score_scalar <- pw * (y - mu) * dmu / vv
  scores <- X * score_scalar

  bread_weight <- pw * dmu^2 / vv
  Xw <- X * sqrt(bread_weight)
  bread <- crossprod(Xw)

  inv_bread <- tryCatch(
    solve(bread),
    error = function(e) NULL
  )

  if (is.null(inv_bread)) {
    stop("The GLM bread matrix is singular; joint covariance is unavailable.",
         call. = FALSE)
  }

  meat <- crossprod(scores)
  vcov_hc0 <- inv_bread %*% meat %*% t(inv_bread)

  dimnames(vcov_hc0) <- list(colnames(X), colnames(X))
  dimnames(inv_bread) <- list(colnames(X), colnames(X))
  colnames(scores) <- colnames(X)

  list(
    frame = mf,
    response = y,
    X = X,
    scores = scores,
    inv_bread = inv_bread,
    vcov_hc0 = vcov_hc0
  )
}

# Internal: joint HC0 covariance of two GLM coefficient vectors.
.baco_joint_glm_covariance <- function(composite_model, critical_model) {
  pc <- .baco_glm_sandwich_parts(composite_model)
  pd <- .baco_glm_sandwich_parts(critical_model)

  rn_c <- rownames(pc$X)
  rn_d <- rownames(pd$X)

  if (!identical(rn_c, rn_d)) {
    stop(
      paste0(
        "The two models must use exactly the same observations, in the same ",
        "order, for joint covariance estimation."
      ),
      call. = FALSE
    )
  }

  if (!identical(colnames(pc$X), colnames(pd$X))) {
    stop(
      paste0(
        "The current implementation requires the same adjustment specification ",
        "in the two models (identical model-matrix columns)."
      ),
      call. = FALSE
    )
  }

  cross_meat <- crossprod(pc$scores, pd$scores)
  cross_vcov <- pc$inv_bread %*% cross_meat %*% t(pd$inv_bread)

  dimnames(cross_vcov) <- list(colnames(pc$X), colnames(pd$X))

  list(
    vcov_composite = pc$vcov_hc0,
    vcov_critical = pd$vcov_hc0,
    cross_vcov = cross_vcov,
    n = nrow(pc$X),
    response_composite = pc$response,
    response_critical = pd$response,
    row_names = rn_c
  )
}

# Internal: use a common row bootstrap and refit both GLMs.
.baco_bootstrap_models <- function(composite_model,
                                   critical_model,
                                   exposure_composite,
                                   exposure_critical,
                                   data,
                                   R,
                                   seed = NULL) {

  if (is.null(data) || !is.data.frame(data)) {
    stop("`data` must be supplied as a data.frame for bootstrap covariance.",
         call. = FALSE)
  }

  if (!is.null(composite_model$call$subset) ||
      !is.null(critical_model$call$subset)) {
    stop(
      paste0(
        "Bootstrap covariance currently does not support models fitted with ",
        "a `subset=` argument. Fit the models on the analysis dataset directly."
      ),
      call. = FALSE
    )
  }

  mf_c <- stats::model.frame(composite_model)
  mf_d <- stats::model.frame(critical_model)

  if (!identical(rownames(mf_c), rownames(mf_d))) {
    stop("The two fitted models must use the same observations.", call. = FALSE)
  }

  used_rows <- rownames(mf_c)
  data_rows <- rownames(data)
  idx <- match(used_rows, data_rows)

  if (anyNA(idx)) {
    stop(
      paste0(
        "Could not map all fitted-model rows back to `data`. Preserve the ",
        "original row names when fitting the models."
      ),
      call. = FALSE
    )
  }

  analysis_data <- data[idx, , drop = FALSE]
  n <- nrow(analysis_data)

  if (length(R) != 1L || !is.finite(R) || R < 100 || abs(R - round(R)) > 1e-8) {
    stop("`R` must be an integer of at least 100 bootstrap replicates.",
         call. = FALSE)
  }
  R <- as.integer(R)

  if (!is.null(seed)) set.seed(seed)

  boot_beta_c <- rep(NA_real_, R)
  boot_beta_d <- rep(NA_real_, R)

  for (b in seq_len(R)) {
    ii <- sample.int(n, size = n, replace = TRUE)
    db <- analysis_data[ii, , drop = FALSE]

    fit_c <- try(
      stats::update(composite_model, data = db),
      silent = TRUE
    )
    fit_d <- try(
      stats::update(critical_model, data = db),
      silent = TRUE
    )

    if (inherits(fit_c, "try-error") || inherits(fit_d, "try-error")) {
      next
    }

    cc <- stats::coef(fit_c)
    cd <- stats::coef(fit_d)

    if (!exposure_composite %in% names(cc) ||
        !exposure_critical %in% names(cd)) {
      next
    }

    bc <- unname(cc[exposure_composite])
    bd <- unname(cd[exposure_critical])

    if (is.finite(bc) && is.finite(bd)) {
      boot_beta_c[b] <- bc
      boot_beta_d[b] <- bd
    }
  }

  ok <- is.finite(boot_beta_c) & is.finite(boot_beta_d)
  n_success <- sum(ok)

  if (n_success < max(50, ceiling(0.50 * R))) {
    stop(
      paste0(
        "Too few successful joint bootstrap refits (", n_success, "/", R, ")."
      ),
      call. = FALSE
    )
  }

  mat <- cbind(
    beta_composite = boot_beta_c[ok],
    beta_critical = boot_beta_d[ok]
  )
  S <- stats::cov(mat)

  list(
    covariance = unname(S[1, 2]),
    se_composite = sqrt(unname(S[1, 1])),
    se_critical = sqrt(unname(S[2, 2])),
    correlation = stats::cor(mat[, 1], mat[, 2]),
    n_success = n_success,
    R = R,
    success_fraction = n_success / R,
    replicates = mat
  )
}

# Internal: infer a user-friendly effect label.
.baco_detect_effect_measure <- function(model_c, model_d) {
  yc <- stats::model.response(stats::model.frame(model_c))
  yd <- stats::model.response(stats::model.frame(model_d))

  binary <- all(yc %in% c(0, 1)) && all(yd %in% c(0, 1))

  if (binary) return("RR")

  fam_c <- model_c$family$family
  fam_d <- model_d$family$family
  if (grepl("poisson", fam_c, ignore.case = TRUE) &&
      grepl("poisson", fam_d, ignore.case = TRUE)) {
    return("IRR")
  }

  "exp(coef)"
}

#' BACO analysis directly from two fitted regression models
#'
#' Extracts the exposure coefficients from two fitted GLMs, estimates their
#' joint sampling covariance, and computes model-based BACO, Fieller inference,
#' Delta, and the joint test of no composite-critical effect discrepancy.
#'
#' @param composite_model Fitted `glm` for the composite endpoint.
#' @param critical_model Fitted `glm` for the critical endpoint.
#' @param exposure Exposure/coefficient name. An exact coefficient name is
#'   preferred. A unique coefficient starting with this string is also accepted.
#' @param covariance_method One of `"joint_sandwich"`, `"bootstrap"`, or
#'   `"both"`. The joint sandwich is the primary analytical method.
#' @param data Original analysis data; required for bootstrap.
#' @param R Number of joint bootstrap replicates.
#' @param seed Optional bootstrap seed.
#' @param conf_level Confidence level.
#' @param effect_measure Optional display label. If `NULL`, the function attempts
#'   to identify `"RR"` for binary log-link GLMs and `"IRR"` for count Poisson
#'   models.
#' @param interpret Logical; if `TRUE`, print a plain-language interpretation
#'   of BACO, Delta, Fieller inference, and the joint discrepancy test.
#'
#' @details
#' Version 0.0.3 intentionally supports GLMs with a log link. The two models
#' must use the same observations and the same adjustment specification.
#'
#' The primary joint covariance is estimated from stacked individual score
#' contributions with an HC0 sandwich covariance. The optional bootstrap uses
#' the same resampled participants for both endpoint models, preserving their
#' dependence.
#'
#' @return An object of class `baco_models_result`.
#' @export
baco_models <- function(composite_model,
                        critical_model,
                        exposure,
                        covariance_method = c("joint_sandwich", "bootstrap", "both"),
                        data = NULL,
                        R = 1000,
                        seed = NULL,
                        conf_level = 0.95,
                        effect_measure = NULL,
                        interpret = FALSE) {

  covariance_method <- match.arg(covariance_method)

  if (!inherits(composite_model, "glm") || !inherits(critical_model, "glm")) {
    stop("The current `baco_models()` implementation supports `glm` objects only.",
         call. = FALSE)
  }

  if (!identical(composite_model$family$link, "log") ||
      !identical(critical_model$family$link, "log")) {
    stop(
      paste0(
        "Both models must currently use a log link so that their coefficients ",
        "represent comparable log-relative effects."
      ),
      call. = FALSE
    )
  }

  coef_c_name <- .baco_resolve_coef(composite_model, exposure)
  coef_d_name <- .baco_resolve_coef(critical_model, exposure)

  beta_c <- unname(stats::coef(composite_model)[coef_c_name])
  beta_d <- unname(stats::coef(critical_model)[coef_d_name])

  if (!is.finite(beta_c) || !is.finite(beta_d)) {
    stop("The requested exposure coefficient was not finite in both models.",
         call. = FALSE)
  }

  joint <- .baco_joint_glm_covariance(composite_model, critical_model)

  # Verify the nested binary endpoint structure whenever this is meaningful.
  yc <- joint$response_composite
  yd <- joint$response_critical

  binary_nested_check <- all(yc %in% c(0, 1)) && all(yd %in% c(0, 1))
  nested_verified <- if (binary_nested_check) {
    all(yd <= yc)
  } else {
    NA
  }

  if (isFALSE(nested_verified)) {
    stop(
      "The critical endpoint is not nested within the composite endpoint (critical <= composite failed).",
      call. = FALSE
    )
  }

  if (is.null(effect_measure)) {
    effect_measure <- .baco_detect_effect_measure(
      composite_model,
      critical_model
    )
  }

  i_c <- match(coef_c_name, rownames(joint$vcov_composite))
  i_d <- match(coef_d_name, rownames(joint$vcov_critical))

  se_c_joint <- sqrt(joint$vcov_composite[i_c, i_c])
  se_d_joint <- sqrt(joint$vcov_critical[i_d, i_d])
  cov_joint <- joint$cross_vcov[i_c, i_d]

  analytic_result <- baco_coef(
    beta_composite = beta_c,
    beta_critical = beta_d,
    se_composite = se_c_joint,
    se_critical = se_d_joint,
    covariance = cov_joint,
    conf_level = conf_level,
    effect_measure = effect_measure,
    interpret = interpret
  )

  bootstrap_info <- NULL
  bootstrap_result <- NULL

  if (covariance_method %in% c("bootstrap", "both")) {
    bootstrap_info <- .baco_bootstrap_models(
      composite_model = composite_model,
      critical_model = critical_model,
      exposure_composite = coef_c_name,
      exposure_critical = coef_d_name,
      data = data,
      R = R,
      seed = seed
    )

    bootstrap_result <- baco_coef(
      beta_composite = beta_c,
      beta_critical = beta_d,
      se_composite = bootstrap_info$se_composite,
      se_critical = bootstrap_info$se_critical,
      covariance = bootstrap_info$covariance,
      conf_level = conf_level,
      effect_measure = effect_measure
    )
  }

  primary_result <- if (identical(covariance_method, "bootstrap")) {
    bootstrap_result
  } else {
    analytic_result
  }

  structure(
    list(
      call = match.call(),
      exposure_requested = exposure,
      coefficient_composite = coef_c_name,
      coefficient_critical = coef_d_name,
      covariance_method = covariance_method,
      n_observations = joint$n,
      nested_verified = nested_verified,
      effect_measure = effect_measure,
      primary = primary_result,
      joint_sandwich = analytic_result,
      bootstrap = bootstrap_result,
      bootstrap_info = bootstrap_info,
      interpret_enabled = isTRUE(interpret)
    ),
    class = "baco_models_result"
  )
}
