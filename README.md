# seurassist

The `seurassist` package provides helper functions for the analysis and
visualization of single-cell RNA sequencing data using Seurat.

## Installation

```r
devtools::install_github("fmbetul/seurassist")
```

## Usage

```r
library(seurassist)
```

Use `?function_name` for full documentation and examples.

### Plotting

| Function | Description |
|---|---|
| `generate_scaled_FPs(seurat_objects, genes, filename, separate)` | Feature plots for genes across several Seurat objects, with the colour scale shared per gene (max expression across objects). Saves PDFs. |
| `generate_AllFPs(SeuFile, gene_list, assay, reduction)` | Feature plots for every gene in a named gene list; saves a PDF. |
| `generate_allVlns(SeuFile, gene_list, assay)` | Violin plots across clusters for every gene in a named gene list; saves a PDF. |
| `make_DotPlot(seurat_obj, gene_list, group_by, assay, ...)` | Grouped `DotPlot` with genes split into named panels; drops genes missing from the assay and de-duplicates across panels. Returns `list(plot, groups_kept)`. |
| `plot_cells_from_sample(sample, ClusterName, integrated, cluster_col)` | Highlights cells of a cluster from a per-sample object on the integrated UMAP. |
| `plot_cells_from_integrated(integrated, ClusterName, sample, cluster_col)` | Highlights cells of a cluster (defined in the integrated object) on a per-sample UMAP. |

### Summaries and clustering helpers

| Function | Description |
|---|---|
| `count_percent_by(object, metadata, colname, sort)` | Cell counts and percentages per value of a metadata column. |
| `choose_dims(stdev, var_threshold, min_contrib, ...)` | Data-driven number of PCA/LSI dimensions from component standard deviations (replaces a fixed `1:30`). |
| `sort_clusters(obj, col)` | Re-levels a cluster column so numeric labels sort numerically (1, 2, ..., 10). |
| `sort_all_clusters(obj, pattern)` | Applies `sort_clusters()` to every metadata column matching a pattern. |
| `save_qc_file(SeuFile, version, qc.plot, Heatmap, csv, ...)` | Saves a UMAP and QC plots to a PDF; optionally a marker heatmap and marker CSVs (`csv = TRUE` requires `Heatmap = TRUE`). |

### QC thresholding

A two-tier approach: tier 1 flags whole bad groups (e.g. clusters); tier 2 sets per-group MAD thresholds for individual cells.

| Function | Description |
|---|---|
| `validate_qc_metrics(metrics, meta)` | Validates a QC metric configuration (a named list of per-column settings such as `direction` and `scale`). |
| `validate_qc_remove_statuses(statuses)` | Validates the group statuses chosen for removal, catching misspellings. |
| `flag_bad_qc_groups(x, group_col, metrics, ...)` | Tier 1: flags groups whose per-metric median is an outlier across groups. |
| `compute_qc_thresholds(x, split_by, metrics, passing_groups, ...)` | Tier 2: per-group, per-metric MAD thresholds for individual cells. |
| `qc_fail_reasons(x, modality, remove_doublets, doublet_col)` | One verdict per cell as a reason string (`""` = passes). |
| `apply_qc_thresholds(seurat_obj, modality, ...)` | Adds `QC_fail_reason` to a Seurat object, subsets to passing cells, and returns summary tables. |

### QC plots

| Function | Description |
|---|---|
| `qc_violin_panel(meta, cluster_col, metric, thresholds, flags, ...)` | Violin of one QC metric per cluster, coloured by tier-1 status with tier-2 thresholds overlaid. |
| `qc_violin_stack(meta, cluster_col, metrics, ...)` | Stacked QC violins, one row per metric. |
| `qc_status_umap(embedding, status, ...)` | UMAP coloured by per-cell QC status. |
| `qc_metric_label(m)` | Display label for a QC metric. |
| `theme_qc(base_size)` | Compact ggplot2 theme for QC plots. |

---

## Dependencies

- Seurat
- dplyr
- ggplot2
- ggpubr
- magrittr
- grDevices
- utils
- scales
- stats
- methods
- rlang

## License

GPL-3

## Issues and Contributions

Please [open an issue](https://github.com/fmbetul/seurassist/issues) on GitHub
for bugs or feature requests.

## Contact

F. M. Betul Erol, MD
