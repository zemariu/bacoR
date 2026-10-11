# bacoR 0.1.0

## First public release

* First public release of `bacoR` on GitHub.
* Added public repository metadata and installation instructions.
* Updated package and citation metadata for the public release.
* Consolidated documentation for aggregate-count, adjusted-effect, fitted-model,
  Fieller, exact-testing, weak-denominator, and zero-event workflows.
* Added links to the companion reproducibility repository for the methodological
  article.
* No statistical algorithms were changed from version 0.0.9.

# bacoR 0.0.9

## R CMD check cleanup

* Added `.Rbuildignore` so development-only `CITATION.cff` and `data-raw/`
  remain in the source repository but are excluded from the built R package.
* Added the standard R package citation file at `inst/CITATION`.
* Validation no longer calls `devtools::build_vignettes()` before
  `devtools::check()`, preventing generated `.R`/`.html` vignette artifacts
  from polluting the source `vignettes/` directory.
* The source package contains only the vignette source `.Rmd`; `R CMD build`
  is allowed to build the vignette normally.
* `interpret = FALSE` remains the default.
* No statistical algorithms were changed from version 0.0.8.

# bacoR 0.0.8

## Documentation and development cleanup

* All 52 unit checks from version 0.0.7 passed.
* `NAMESPACE` is now explicitly managed by roxygen2.
* Removed first-pass unresolved documentation links from the package overview.
* The validation script now clears previously sourced `baco_*` functions from
  the global environment before `load_all()`, preventing misleading masking
  warnings during development.
* `interpret = FALSE` remains the default in all public analytical interfaces.
* Removed stale version-specific wording from `baco_models()` error messages.

# bacoR 0.0.7

## Packaging fix

* Fixed the DCF structure of `DESCRIPTION`; blank lines are no longer present.
* Added `knitr` and `rmarkdown` correctly under `Suggests`.
* Confirmed `interpret = FALSE` remains the default in all public analysis
  functions.
* Retained the introductory vignette and package-level documentation from
  version 0.0.6.

# bacoR 0.0.6

## Documentation milestone

* Added a package-level overview and a complete introductory vignette.
* Confirmed `interpret = FALSE` as the default in all public analytical
  interfaces.
* Documented aggregate-count, adjusted-effect, fitted-model, bootstrap,
  zero-event, and Fieller workflows.
* Added recommended reporting guidance.

# bacoR 0.0.5

* Added optional `interpret = TRUE`.
* BACO point-estimate interpretation follows the original Diaz-Quijano
  classification.
* Added interpretation for Delta/QRR, discrepancy testing, and Fieller
  diagnostics.

# bacoR 0.0.4

* Added the simulated `baco_trial` dataset.
* Simplified standard output and Fieller reporting.
