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

---

### 1. Generate Scaled Feature Plots

Generates feature plots for a list of genes across multiple Seurat objects,
scaling the colour axis by the maximum expression value of each gene across
all objects. Saves the result as PDF files.

```r
generate_scaled_FPs(seurat_objects, genes, filename, separate = TRUE)
```

| Argument | Description |
|---|---|
| `seurat_objects` | A named list of Seurat objects. |
| `genes` | A character vector of gene names. |
| `filename` | Output filename (without extension). Defaults to today's date. |
| `separate` | If `TRUE`, saves a separate PDF per Seurat object. Default is `TRUE`. |

---

### 2. Generate Feature Plots for Given Genes

Generates feature plots for each gene in `gene_list` and saves them as a PDF.

```r
generate_AllFPs(SeuFile, gene_list, assay = "RNA", reduction = "umap")
```

| Argument | Description |
|---|---|
| `SeuFile` | A Seurat object. |
| `gene_list` | A named list of character vectors of gene names. |
| `assay` | Assay to use. Default is `"RNA"`. |
| `reduction` | Dimensional reduction to use. Default is `"umap"`. |

---

### 3. Generate Violin Plots for Given Genes Across Clusters

Generates violin plots for each gene in `gene_list` across clusters and saves
them as a PDF.

```r
generate_allVlns(SeuFile, gene_list, assay = "RNA")
```

| Argument | Description |
|---|---|
| `SeuFile` | A Seurat object. |
| `gene_list` | A named list of character vectors of gene names. |
| `assay` | Assay to use. Default is `"RNA"`. |

---

### 4. Save QC File

Saves a UMAP plot and quality control (QC) plots to a PDF. Optionally generates
a marker gene heatmap and exports CSV files.

> **Note:** `csv = TRUE` requires `Heatmap = TRUE`, since marker genes must be
> computed before they can be exported.

```r
save_qc_file(SeuFile, version = "v1.0.0", qc.plot = TRUE, Heatmap = FALSE,
             csv = FALSE, pdf.height = 10, pdf.width = 12)
```

| Argument | Description |
|---|---|
| `SeuFile` | A Seurat object. |
| `version` | Version string used in the output filename. Default is `"v1.0.0"`. |
| `qc.plot` | If `TRUE`, generates per-cluster QC violin plots. Default is `TRUE`. |
| `Heatmap` | If `TRUE`, generates a top-marker heatmap. Default is `FALSE`. |
| `csv` | If `TRUE` and `Heatmap = TRUE`, exports marker CSVs. Default is `FALSE`. |
| `pdf.height` | PDF height in inches. Default is `10`. |
| `pdf.width` | PDF width in inches. Default is `12`. |

---

### 5. Count and Percentage Summary of Metadata

Returns a data frame with cell counts and percentages for each unique value in
a metadata column. Useful for summarising cell type or sample compositions.

```r
count_percent_by(object, metadata, colname = "CellType", sort = TRUE)
```

| Argument | Description |
|---|---|
| `object` | A Seurat object. |
| `metadata` | Name of the metadata column to summarise (e.g. `"celltype"`, `"sample"`). |
| `colname` | Column name for the category in the output data frame. Default is `"CellType"`. |
| `sort` | If `TRUE`, sorts by count in decreasing order. Default is `TRUE`. |

**Example output:**

```
    CellType Count Percent
1    T cells  1200    35.2
2   B cells    850    24.9
3  Monocytes    600    17.6
...
```

---

### 6. Plot Cells of a Specified Sample and Cluster on the Integrated UMAP

Highlights cells from a specified cluster in a sample on the UMAP of an
integrated Seurat object.

```r
plot_cells_from_sample(sample, ClusterName, integrated, cluster_col = "celltype")
```

| Argument | Description |
|---|---|
| `sample` | A per-sample Seurat object. |
| `ClusterName` | The cell type or cluster name to highlight. |
| `integrated` | The integrated Seurat object. |
| `cluster_col` | Metadata column holding cluster labels. Default is `"celltype"`. |

---

### 7. Plot Specific Cells from an Integrated Object on a Per-Sample UMAP

Highlights cells from a specified cluster (as defined in the integrated object)
on the UMAP of an individual sample.

```r
plot_cells_from_integrated(integrated, ClusterName, sample, cluster_col = "celltype")
```

| Argument | Description |
|---|---|
| `integrated` | The integrated Seurat object. |
| `ClusterName` | The cell type or cluster name to highlight. |
| `sample` | A per-sample Seurat object. |
| `cluster_col` | Metadata column holding cluster labels. Default is `"celltype"`. |

---

## Dependencies

- Seurat
- dplyr
- ggplot2
- ggpubr
- magrittr
- grDevices
- utils
- rlang

## License

GPL-3

## Issues and Contributions

Please [open an issue](https://github.com/fmbetul/seurassist/issues) on GitHub
for bugs or feature requests.

## Contact

F. M. Betul Erol — fmab_erol@hotmail.com
