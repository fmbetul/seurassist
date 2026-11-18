# QC plotting helpers for the MAD-threshold workflow (qc_thresholds.R).
# Adapted from the 28-novo-in-vivo pipeline.

# Level-1 group status (violins)
QC_STATUS_COLORS <- c(pass             = "#2e7d32",
                      review           = "#e65100",
                      review_small_n   = "#6a1b9a",
                      low_qc           = "#c62828",
                      possible_doublet = "#00838f")

# Per-cell QC verdict (UMAP)
QC_CELL_COLORS <- c(pass = "grey78", fail_RNA = "#d32f2f",
                    fail_ATAC = "#ffb300", fail_both = "#4a148c",
                    doublet = "#1565c0")

THR_MIN_COLOR <- "#d32f2f"
THR_MAX_COLOR <- "#1565c0"

# Metrics drawn on a log10 axis (counts span orders of magnitude)
QC_LOG10_METRICS <- c("nCount_RNA", "nFeature_RNA",
                      "nCount_ATAC_cCRE", "nFeature_ATAC_cCRE",
                      "atac_fragments",
                      "nCount_CRISPR", "nFeature_CRISPR")

QC_METRIC_LABELS <- c(
  nCount_RNA         = "RNA UMIs",
  nFeature_RNA       = "RNA genes",
  percent.mt         = "MT %",
  nCount_ATAC_cCRE   = "ATAC counts in cCREs",
  nFeature_ATAC_cCRE = "cCREs detected",
  TSS.enrichment     = "TSS enrichment",
  nucleosome_signal  = "Nucleosome signal",
  blacklist_ratio    = "Blacklist ratio",
  frip_ccre          = "FRiP (cCRE)",
  atac_fragments     = "ATAC fragments",
  scDblFinder.score  = "Doublet score",
  nCount_CRISPR      = "Guide UMIs (all features)",
  nFeature_CRISPR    = "Guide features detected")

# `QC_METRIC_LABELS[m]` returns NA (not NULL) for an m absent from the map, so
# `%||%` never sees the NULL it's looking for and the panel title silently
# reads "NA" -- this bit every metric added after the map was last updated.
#' Display label for a QC metric
#'
#' @param m Character vector of metric (column) names.
#' @return Human-readable labels; names not in `QC_METRIC_LABELS` are returned
#'   unchanged.
#' @export
qc_metric_label <- function(m) {
  lbl <- unname(QC_METRIC_LABELS[m])
  ifelse(is.na(lbl), m, lbl)
}

#' Compact ggplot2 theme for QC plots
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @export
theme_qc <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(size = base_size + 1, face = "bold"),
          legend.key.size  = unit(0.8, "lines"))
}

#' Violin of one QC metric per cluster with tier-1 status and tier-2 thresholds
#'
#' Violins are coloured by tier-1 group status, with the tier-2 min (red) and max
#' (blue) threshold lines drawn on the clusters that have them. A red violin says
#' tier 1 wants the whole cluster gone; a threshold line marks the cells tier 2
#' removes from a cluster that stays. A thin boxplot (median/IQR) is overlaid.
#'
#' @param meta Cell metadata.
#' @param cluster_col Cluster column plotted on x.
#' @param metric Metric column plotted on y.
#' @param thresholds [compute_qc_thresholds()] output, or `NULL`.
#' @param flags [flag_bad_qc_groups()] output, or `NULL`.
#' @param log10_y Use a log10 y axis; defaults to metrics in `QC_LOG10_METRICS`.
#' @param show_legend Show the group-status legend.
#' @param legend_direction `"horizontal"` or `"vertical"`.
#' @param legend_text_size `NULL` or a point size for legend text and title.
#' @param title_size `NULL` or a point size for the panel title.
#' @return A ggplot object.
#' @examples
#' set.seed(1)
#' meta <- data.frame(cl = rep(1:4, each = 50), nCount_RNA = rpois(200, 3000))
#' mt <- list(nCount_RNA = list(direction = "low_bad"))
#' qc_violin_panel(meta, "cl", "nCount_RNA",
#'                 thresholds = compute_qc_thresholds(meta, "cl", mt),
#'                 flags = flag_bad_qc_groups(meta, "cl", mt),
#'                 show_legend = TRUE)
#' @export
qc_violin_panel <- function(meta, cluster_col, metric,
                            thresholds = NULL, flags = NULL,
                            log10_y = metric %in% QC_LOG10_METRICS,
                            show_legend = FALSE,
                            legend_direction = "horizontal",
                            legend_text_size = NULL,
                            title_size = NULL) {

  lv <- sorted_levels(meta[[cluster_col]])
  df <- data.frame(cluster = factor(as.character(meta[[cluster_col]]), levels = lv),
                   value   = meta[[metric]],
                   stringsAsFactors = FALSE)
  df <- df[is.finite(df$value), , drop = FALSE]

  if (!is.null(flags)) {
    st <- stats::setNames(flags$status, as.character(flags$group))
    df$status <- factor(unname(st[as.character(df$cluster)]),
                        levels = names(QC_STATUS_COLORS))
  } else {
    df$status <- factor("pass", levels = names(QC_STATUS_COLORS))
  }

  p <- ggplot(df, aes(x = .data$cluster, y = .data$value, fill = .data$status)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey30",
               alpha = 0.72, trim = TRUE, na.rm = TRUE, show.legend = FALSE) +
    geom_boxplot(width = 0.08, fill = "white", colour = "grey20",
                coef = 0, outlier.shape = NA, linewidth = 0.35,
                na.rm = TRUE, show.legend = FALSE) +
    scale_fill_manual(values = QC_STATUS_COLORS, drop = FALSE,
                      name = "group status",
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    labs(x = cluster_col, y = qc_metric_label(metric),
         title = qc_metric_label(metric)) +
    theme_qc() +
    theme(legend.position  = if (show_legend) "top" else "none",
          legend.direction = legend_direction,
          legend.key.size  = unit(8, "pt"))

  if (!is.null(legend_text_size))
    p <- p + theme(legend.text  = element_text(size = legend_text_size),
                   legend.title = element_text(size = legend_text_size))
  if (!is.null(title_size))
    p <- p + theme(plot.title = element_text(size = title_size, face = "bold"))

  # geom_violin/geom_boxplot run stats (stat_ydensity/stat_boxplot) that drop
  # a status level with zero cells in this sample *before* legend keys are
  # built, so scale_fill_manual(drop = FALSE) alone still leaves an absent
  # status without a colour box. A single invisible geom_point layer
  # (stat = "identity", so nothing gets dropped) carries every level instead,
  # sized/alpha'd to zero so it has no visual effect on the panel itself.
  if (show_legend && !is.null(flags)) {
    legend_df <- data.frame(
      cluster = factor(lv[1], levels = lv),
      value   = min(df$value, na.rm = TRUE),
      status  = factor(names(QC_STATUS_COLORS), levels = names(QC_STATUS_COLORS)),
      stringsAsFactors = FALSE)
    p <- p + geom_point(data = legend_df, aes(fill = .data$status),
                        alpha = 0, size = 0, stroke = 0, na.rm = TRUE,
                        show.legend = TRUE, key_glyph = "rect")
  }

  if (log10_y)
    p <- p + scale_y_log10(labels = scales::label_comma()) +
      annotation_logticks(sides = "l", linewidth = 0.25)

  if (!is.null(thresholds)) {
    t <- thresholds[thresholds$metric == metric, , drop = FALSE]
    if (nrow(t)) {
      t$x <- match(as.character(t$group), lv)
      t   <- t[!is.na(t$x), , drop = FALSE]
      seg <- function(col, colour) {
        d <- t[is.finite(t[[col]]) & (!log10_y | t[[col]] > 0), , drop = FALSE]
        if (!nrow(d)) return(NULL)
        geom_segment(data = d, inherit.aes = FALSE,
                     aes(x = .data$x - 0.45, xend = .data$x + 0.45,
                         y = .data[[col]], yend = .data[[col]]),
                     colour = colour, linewidth = 0.5)
      }
      p <- p + seg("min_thr", THR_MIN_COLOR) + seg("max_thr", THR_MAX_COLOR)
    }
  }
  p
}


#' Stacked QC violins, one row per metric
#'
#' All clusters of one clustering side by side on a shared x axis, one facet row
#' per metric. Meant to be called once per page, not assembled like
#' [qc_violin_panel()]'s small multiples.
#'
#' @param meta Cell metadata.
#' @param cluster_col Cluster column plotted on x.
#' @param metrics Character vector of metric ids, top to bottom. An id ending in
#'   `"_log2"` plots `log2()` of the column named by stripping that suffix.
#' @param title Plot title.
#' @param x_angle Rotation of x-axis tick labels (45 suits long labels).
#' @param split_by Optional metadata column; when set, violins are dodged within
#'   each cluster and coloured by it.
#' @return A ggplot object, or `NULL` if no requested metric has finite values.
#' @examples
#' set.seed(1)
#' meta <- data.frame(cl = rep(1:4, each = 50), nCount_RNA = rpois(200, 3000),
#'                    percent.mt = runif(200, 0, 5))
#' qc_violin_stack(meta, "cl", c("nCount_RNA", "nCount_RNA_log2", "percent.mt"))
#' @export
qc_violin_stack <- function(meta, cluster_col, metrics, title = NULL, x_angle = 0,
                            split_by = NULL) {

  lv <- sorted_levels(meta[[cluster_col]])

  row_label <- function(m)
    if (endsWith(m, "_log2")) paste0("log2 ", qc_metric_label(sub("_log2$", "", m)))
    else qc_metric_label(m)

  rows <- lapply(metrics, function(m) {
    base <- if (endsWith(m, "_log2")) sub("_log2$", "", m) else m
    if (!base %in% colnames(meta)) return(NULL)
    val <- meta[[base]]
    if (endsWith(m, "_log2")) val <- log2(val)
    out <- data.frame(cluster = factor(as.character(meta[[cluster_col]]), levels = lv),
                      metric  = m, value = val, stringsAsFactors = FALSE)
    if (!is.null(split_by)) out[[split_by]] <- meta[[split_by]]
    out
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  df <- df[is.finite(df$value), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  present <- intersect(metrics, unique(df$metric))
  df$metric <- factor(df$metric, levels = present,
                      labels = vapply(present, row_label, character(1)))

  p <- if (is.null(split_by)) {
    ggplot(df, aes(x = .data$cluster, y = .data$value, fill = .data$metric)) +
      geom_violin(scale = "width", linewidth = 0.2, colour = "grey20",
                 trim = TRUE, na.rm = TRUE, show.legend = FALSE)
  } else {
    ggplot(df, aes(x = .data$cluster, y = .data$value, fill = .data[[split_by]])) +
      geom_violin(scale = "width", linewidth = 0.2, colour = "grey20",
                 trim = TRUE, na.rm = TRUE, position = position_dodge(width = 0.9))
  }

  p +
    facet_grid(rows = vars(.data$metric), scales = "free_y") +
    labs(x = cluster_col, y = "value", title = title) +
    theme_qc() +
    theme(strip.text.y      = element_text(angle = 0, hjust = 0),
          panel.spacing.y   = unit(2, "pt"),
          axis.text.x       = element_text(angle = x_angle,
                                           hjust = if (x_angle == 0) 0.5 else 1))
}


#' UMAP scatter coloured by per-cell QC status
#'
#' Built directly from an embedding matrix rather than a Seurat object, so it
#' works on a lightweight metadata + embedding table where [Seurat::DimPlot()]
#' is not usable.
#'
#' @param embedding Numeric matrix/data.frame, cells x 2, e.g.
#'   `Seurat::Embeddings(obj, "umap")`.
#' @param status Factor/character of per-cell QC status, same length as
#'   `nrow(embedding)`, with levels matching `names(QC_CELL_COLORS)`.
#' @param title Plot title.
#' @param pt_size Point size.
#' @param order Draw non-`"pass"` cells last (on top of the background).
#' @param legend_text_size,legend_title_size Legend label and title sizes.
#' @param legend_key_point_size Size of the glyph in each legend key.
#' @return A ggplot object.
#' @examples
#' set.seed(1)
#' emb <- matrix(rnorm(400), 200, 2, dimnames = list(NULL, c("UMAP_1", "UMAP_2")))
#' status <- factor(sample(c("pass", "fail_RNA"), 200, TRUE, c(.9, .1)),
#'                  levels = names(QC_CELL_COLORS))
#' qc_status_umap(emb, status, title = "umap")
#' @export
qc_status_umap <- function(embedding, status, title = NULL, pt_size = 0.2,
                           order = TRUE, legend_text_size = 7.7,
                           legend_title_size = 8.4, legend_key_point_size = 2.8) {
  df <- data.frame(x = embedding[, 1], y = embedding[, 2],
                   status = status, stringsAsFactors = FALSE)
  if (order) df <- df[order(df$status == "pass", decreasing = TRUE), , drop = FALSE]
  ggplot(df, aes(x = .data$x, y = .data$y, colour = .data$status)) +
    geom_point(size = pt_size, alpha = 0.5, na.rm = TRUE) +
    scale_colour_manual(values = QC_CELL_COLORS, drop = FALSE, name = "QC_preview",
                        guide = guide_legend(override.aes = list(
                          size = legend_key_point_size, alpha = 1))) +
    labs(x = colnames(embedding)[1], y = colnames(embedding)[2], title = title) +
    theme_qc() +
    theme(legend.text     = element_text(size = legend_text_size),
          legend.title    = element_text(size = legend_title_size),
          legend.key.size = unit(1.05, "lines"))
}



#' QC constants: status levels, colours and metric labels
#'
#' Constants used by the QC threshold and plotting functions.
#'
#' \describe{
#'   \item{`QC_GROUP_STATUSES`}{Group statuses [flag_bad_qc_groups()] can assign.}
#'   \item{`QC_REMOVE_STATUSES_DEFAULT`}{Statuses that remove a whole group by default.}
#'   \item{`QC_STATUS_COLORS`}{Fill colours for group status.}
#'   \item{`QC_CELL_COLORS`}{Colours for per-cell QC status in [qc_status_umap()].}
#'   \item{`THR_MIN_COLOR`, `THR_MAX_COLOR`}{Colours of min/max threshold lines.}
#'   \item{`QC_LOG10_METRICS`}{Metrics drawn on a log10 axis by [qc_violin_panel()].}
#'   \item{`QC_METRIC_LABELS`}{Named vector of display labels, see [qc_metric_label()].}
#' }
#'
#' @name qc_constants
#' @aliases QC_GROUP_STATUSES QC_REMOVE_STATUSES_DEFAULT QC_STATUS_COLORS
#'   QC_CELL_COLORS THR_MIN_COLOR THR_MAX_COLOR QC_LOG10_METRICS QC_METRIC_LABELS
NULL
