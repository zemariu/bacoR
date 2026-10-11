# bacoR: Robust inference for Bias Attributable to Composite Outcomes

[![R-CMD-check](https://github.com/zemariu/bacoR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zemariu/bacoR/actions/workflows/R-CMD-check.yaml)
[![GitHub release](https://img.shields.io/github/v/release/zemariu/bacoR)](https://github.com/zemariu/bacoR/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md)

`bacoR` is an R package for estimation and inference for the **Bias Attributable
to Composite Outcome (BACO)** framework. It provides tools for ratio-based BACO
estimation, direct composite-critical effect discrepancy measures, Fieller
confidence sets, exact conditional inference for nested aggregate binary
outcomes, and model-based analyses using joint sampling covariance.

Version **0.1.0** is the first public release.

## Main functionality

The package provides:

- conventional BACO point estimation when mathematically defined;
- Fieller confidence sets for the BACO ratio;
- the direct log-effect discrepancy
  `Delta = log(RR_composite) - log(RR_critical)`;
- `QRR = (D_t/C_t)/(D_c/C_c) = exp(-Delta)`;
- the probability-ordered two-sided Fisher exact test of no
  composite-critical effect discrepancy for nested aggregate binary outcomes;
- optional probability-ordered Fisher mid-P sensitivity analysis;
- diagnostics for zero critical-event cells and weak-denominator behavior;
- inference from estimated log effects and their joint covariance;
- model-based inference from fitted log-link models using joint sandwich
  covariance, with optional joint bootstrap support;
- no automatic continuity correction.

The main exported functions are `baco_counts()`, `baco_fieller()`,
`baco_coef()`, `baco_effects()`, and `baco_models()`. The package also includes
the simulated `baco_trial` dataset for examples, teaching, and software testing.

## Installation

Install the public development release from GitHub with:

```r
install.packages("remotes")
remotes::install_github("zemariu/bacoR")
```

Then load the package:

```r
library(bacoR)
```

## Aggregate-count analysis

For a two-group comparison with a binary composite endpoint and a prespecified
critical component nested within that composite, use `baco_counts()`:

```r
x <- baco_counts(
  n_treatment = 1000,
  n_control = 1000,
  composite_treatment = 150,
  composite_control = 200,
  critical_treatment = 30,
  critical_control = 60
)

x
```

The output includes the BACO point estimate when defined, Fieller ratio
inference, `Delta`, QRR, weak-denominator diagnostics, and the exact
composite-critical discrepancy test.

For count-based nested binary outcomes, the probability-ordered two-sided
Fisher exact test is the primary discrepancy test. The probability-ordered
mid-P version is available as a sensitivity analysis by setting `midp = TRUE`.
Zero critical-event cells are retained without automatic continuity correction.

## Inference from adjusted effects

When the composite and critical effects are estimated from the same
participants, their sampling dependence must be accounted for. If the joint
covariance is available, use `baco_coef()`:

```r
x <- baco_coef(
  beta_composite = log(0.80),
  beta_critical = log(0.60),
  se_composite = 0.08,
  se_critical = 0.12,
  covariance = 0.00576,
  effect_measure = "RR"
)

x
```

The covariance cannot generally be reconstructed from the two standard errors
alone.

If adjusted RR, HR, or IRR estimates and the dependence between their log
estimates are available, `baco_effects()` provides a direct interface:

```r
baco_effects(
  effect_composite = 0.80,
  effect_critical = 0.60,
  se_log_composite = 0.08,
  se_log_critical = 0.12,
  correlation_log = 0.60,
  effect_measure = "RR"
)
```

## Inference from fitted models

`baco_models()` can be used with two compatible log-link models fitted to the
same participants and with the same adjustment specification:

```r
fit_comp <- glm(
  composite ~ treatment + age + comorbidity,
  family = poisson(link = "log"),
  data = dat
)

fit_crit <- glm(
  critical ~ treatment + age + comorbidity,
  family = poisson(link = "log"),
  data = dat
)

baco_models(
  composite_model = fit_comp,
  critical_model = fit_crit,
  exposure = "treatment",
  covariance_method = "joint_sandwich",
  effect_measure = "RR"
)
```

The function extracts the exposure coefficients and estimates their dependence
using joint participant-level score contributions. Joint bootstrap support is
also available for model-based sensitivity or covariance validation. Bootstrap
inference was not evaluated as a competing primary inferential procedure in the
methodological study associated with this package.

## Fieller confidence sets and weak denominators

A bounded Fieller set is a finite interval for the BACO ratio. Disjoint
unbounded and all-real sets are retained explicitly rather than being forced
into conventional finite confidence intervals. These non-bounded forms are
important diagnostics of weak ratio identification when the critical-endpoint
log effect is close to zero.

## Simulated example dataset

`baco_trial` is an individual-level simulated trial included for examples,
teaching, vignettes, and software testing. It contains no real or confidential
participant data.

```r
head(baco_trial)
```

## Optional interpretation

The main analysis interfaces accept `interpret = TRUE`. This adds a structured
plain-language interpretation while leaving all numerical results unchanged.
Statistical testing is deliberately separated from descriptive BACO
classification: a non-significant discrepancy test indicates absence of
evidence of a difference, not evidence of equivalence.

## Reproducibility materials

The simulation scripts, empirical reanalysis workflow, standardized data, and
frozen results associated with the methodological article are maintained in a
separate repository:

<https://github.com/zemariu/BACO_article>

## Associated methodological study

The package was developed alongside the manuscript:

> *Wald inference for the BACO ratio is poorly calibrated when the critical
> effect is near the null: a simulation and empirical study.*

The reproducibility materials associated with the manuscript are archived on
Zenodo (DOI: [10.5281/zenodo.23290575](https://doi.org/10.5281/zenodo.23290575)).

The `preferred-citation` metadata in `CITATION.cff` will be updated with the
article DOI and publication details after publication.

## License

`bacoR` is distributed under the GNU General Public License version 3 or later
(GPL-3.0-or-later).

Year: 2026

Copyright Holder: bacoR authors
