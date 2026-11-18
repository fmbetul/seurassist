# seurassist 0.2.0

First release under the name `seurassist` (formerly `bioinfotools`).

## New
* Two-tier QC thresholding: `flag_bad_qc_groups()`, `compute_qc_thresholds()`,
  `qc_fail_reasons()`, `apply_qc_thresholds()`, with input validators
  `validate_qc_metrics()` and `validate_qc_remove_statuses()`.
* QC plots: `qc_violin_panel()`, `qc_violin_stack()`, `qc_status_umap()`,
  `qc_metric_label()`, `theme_qc()`.
* `choose_dims()` for data-driven PCA/LSI dimension selection.
* `sort_clusters()` / `sort_all_clusters()` for numeric-aware cluster ordering.
* `make_DotPlot()` for grouped dot plots with named gene panels.
* `count_percent_by()` for metadata count/percentage summaries.
* A testthat test suite.

## Fixes
* Bug fixes in `generate_scaled_FPs()` and `make_DotPlot()`.
* Removed non-ASCII characters from code and silenced undefined-variable notes
  reported by `R CMD check`.

## Other
* README rewritten to cover all exported functions.
