# ==============================================================================
# Build the simulated baco_trial dataset
#
# Development-only script for the bacoR source repository.
# The data-raw/ directory is excluded from the built package through
# .Rbuildignore.
#
# Running this script from the package root recreates:
#
#   data/baco_trial.rda
#
# The dataset is simulated and contains no real or confidential participant
# information.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Generate the dataset
# ------------------------------------------------------------------------------

generate_baco_trial <- function(seed = 20260902L, n = 2500L) {

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)

  treatment <- stats::rbinom(n, 1, 0.50)

  age <- stats::rnorm(n, mean = 50, sd = 11)
  age <- pmin(pmax(age, 18), 85)

  sex <- stats::rbinom(n, 1, 0.50)

  comorbidity <- stats::rbinom(
    n,
    1,
    stats::plogis(
      -1.15 +
        0.018 * (age - 50) +
        0.25 * sex
    )
  )

  lp_composite <-
    -1.45 +
    log(0.80) * treatment +
    0.006 * (age - 50) +
    0.16 * sex +
    0.24 * comorbidity

  lp_critical <-
    -2.80 +
    log(0.60) * treatment +
    0.006 * (age - 50) +
    0.16 * sex +
    0.24 * comorbidity

  p_composite <- exp(lp_composite)
  p_critical <- exp(lp_critical)

  if (
    any(p_composite >= 1) ||
    any(p_critical >= 1) ||
    any(p_critical > p_composite)
  ) {
    stop("Internal error while generating `baco_trial`.")
  }

  # A shared uniform variate guarantees critical <= composite.
  u <- stats::runif(n)

  critical <- as.integer(u < p_critical)
  composite <- as.integer(u < p_composite)

  out <- data.frame(
    id = seq_len(n),
    treatment = treatment,
    age = round(age, 1),
    sex = sex,
    comorbidity = comorbidity,
    composite = composite,
    critical = critical
  )

  attr(out, "data_generating_parameters") <- list(
    rr_composite = 0.80,
    rr_critical = 0.60,
    baco = log(0.80) / log(0.60),
    delta = log(0.80) - log(0.60)
  )

  out
}

baco_trial <- generate_baco_trial()


# ------------------------------------------------------------------------------
# 2. Validate the frozen dataset
# ------------------------------------------------------------------------------

expected_names <- c(
  "id",
  "treatment",
  "age",
  "sex",
  "comorbidity",
  "composite",
  "critical"
)

stopifnot(
  is.data.frame(baco_trial),
  nrow(baco_trial) == 2500L,
  identical(names(baco_trial), expected_names),
  identical(baco_trial$id, seq_len(nrow(baco_trial))),
  all(baco_trial$treatment %in% c(0L, 1L)),
  all(baco_trial$sex %in% c(0L, 1L)),
  all(baco_trial$comorbidity %in% c(0L, 1L)),
  all(baco_trial$composite %in% c(0L, 1L)),
  all(baco_trial$critical %in% c(0L, 1L)),
  all(baco_trial$critical <= baco_trial$composite)
)

pars <- attr(baco_trial, "data_generating_parameters")

stopifnot(
  is.list(pars),
  isTRUE(all.equal(pars$rr_composite, 0.80)),
  isTRUE(all.equal(pars$rr_critical, 0.60)),
  isTRUE(all.equal(pars$baco, log(0.80) / log(0.60))),
  isTRUE(all.equal(pars$delta, log(0.80) - log(0.60)))
)


# ------------------------------------------------------------------------------
# 3. Save package data
# ------------------------------------------------------------------------------

dir.create("data", showWarnings = FALSE, recursive = TRUE)

save(
  baco_trial,
  file = file.path("data", "baco_trial.rda"),
  compress = "xz",
  version = 2
)

message(
  "Created data/baco_trial.rda (",
  nrow(baco_trial),
  " rows, ",
  ncol(baco_trial),
  " variables)."
)
