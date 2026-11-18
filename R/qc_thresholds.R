# Two-tier, MAD-based QC thresholding (see header of each function).
# Adapted from the 28-novo-in-vivo pipeline.

# Group statuses tier 1 can assign. Anything outside this set is a typo, and a
# typo in `remove_statuses` would silently retain a cluster meant for removal.
QC_GROUP_STATUSES <- c("pass", "review", "review_small_n",
                       "low_qc", "possible_doublet")

# Statuses that remove the whole group. `review` is included: tier 1 flagged it
# on one metric, and a cluster that is an outlier on one metric is not something
# to keep silently. Override per script if that is too aggressive -- but
# override in one place, and the value is written into the flags CSV.
QC_REMOVE_STATUSES_DEFAULT <- c("low_qc", "possible_doublet",
                                "review", "review_small_n")

#' Validate a set of group statuses to remove
#'
#' A misspelt status matches no group, so the clusters it was meant to remove are
#' silently kept. Call this wherever the removal set is configured.
#'
#' @param statuses Character vector of statuses; must be a subset of
#'   `QC_GROUP_STATUSES` (`pass`, `review`, `review_small_n`, `low_qc`,
#'   `possible_doublet`).
#' @param what Label used in error messages.
#' @return `statuses` (unique), invisibly. Errors on unknown statuses.
#' @examples
#' validate_qc_remove_statuses(c("low_qc", "possible_doublet"))
#' @export
validate_qc_remove_statuses <- function(statuses,
                                        what = "remove_statuses") {
  if (!is.character(statuses))
    stop(what, " must be a character vector", call. = FALSE)
  bad <- setdiff(statuses, QC_GROUP_STATUSES)
  if (length(bad))
    stop(what, " contains status(es) flag_bad_qc_groups() never assigns: ",
         paste(bad, collapse = ", "),
         "\n  known statuses: ", paste(QC_GROUP_STATUSES, collapse = ", "),
         call. = FALSE)
  invisible(unique(statuses))
}


# ── Internal helpers ─────────────────────────────────────────────────────────

# meta.data from a Seurat object, or the data.frame itself. Keeps every
# decision function usable in tests with no Seurat installed.
.qc_meta <- function(x) {
  if (is.data.frame(x)) return(x)
  if (isTRUE(methods::is(x, "Seurat"))) return(x@meta.data)
  stop("expected a data.frame of cell metadata or a Seurat object, got: ",
       paste(class(x), collapse = "/"), call. = FALSE)
}

.to_scaled <- function(x, scale) {
  if (scale != "log2") return(x)
  if (any(x < -1, na.rm = TRUE))
    stop("scale = \"log2\" was requested for values below -1; ",
         "log2(x + 1) is undefined there. Use scale = \"raw\" for signed ",
         "metrics. Offending minimum: ", format(min(x, na.rm = TRUE)),
         call. = FALSE)
  log2(x + 1)
}

.from_scaled <- function(x, scale) if (scale == "log2") 2^x - 1 else x

# Cluster labels sort numerically when they look numeric ("2" before "10"),
# lexically otherwise. Purely cosmetic, but the threshold CSVs are read by a
# human reviewing the threshold tables.
.sort_groups <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  n <- suppressWarnings(as.numeric(x))
  if (!anyNA(n)) x[order(n)] else sort(x)
}

.finite <- function(x) x[is.finite(x)]

# Clamp a threshold to the biological guardrails. A guardrail is a limit the
# algorithm may not exceed, never the threshold itself.
.apply_guardrails <- function(thr, direction, hard_floor, hard_ceiling) {
  if (direction == "low_bad") {
    if (!is.null(hard_floor) && is.finite(hard_floor)) thr <- max(thr, hard_floor)
  } else {
    if (!is.null(hard_ceiling) && is.finite(hard_ceiling)) thr <- min(thr, hard_ceiling)
    if (!is.null(hard_floor) && is.finite(hard_floor)) thr <- max(thr, hard_floor)
  }
  thr
}

# One metric, one group. Returns the bounds plus enough provenance to see which
# rule actually set them.
.compute_one <- function(x, cfg, upper_quantile = NULL, min_cells = 5L) {

  direction <- cfg$direction
  scale     <- cfg$scale %||% "log2"
  n_mad     <- cfg$n_mad %||% 3

  x <- .finite(x)
  if (length(x) < min_cells)
    return(list(min_thr = NA_real_, max_thr = NA_real_, max_thr_mad = NA_real_,
                median_raw = NA_real_, mad_raw = NA_real_,
                cap_applied = FALSE, method = "too_few_cells"))

  xs    <- .to_scaled(x, scale)
  med_s <- stats::median(xs)
  mad_s <- stats::mad(xs)

  cap <- if (is.null(upper_quantile)) NA_real_ else
    as.numeric(stats::quantile(x, probs = upper_quantile, names = FALSE))

  if (direction == "low_bad") {
    min_thr     <- .apply_guardrails(.from_scaled(med_s - n_mad * mad_s, scale),
                                     direction, cfg$hard_floor, cfg$hard_ceiling)
    max_thr_mad <- NA_real_
    max_thr     <- cap                       # only the cap bounds a low_bad metric
    cap_applied <- is.finite(cap)
  } else {
    min_thr     <- NA_real_
    max_thr_mad <- .apply_guardrails(.from_scaled(med_s + n_mad * mad_s, scale),
                                     direction, cfg$hard_floor, cfg$hard_ceiling)
    max_thr     <- if (is.finite(cap)) min(max_thr_mad, cap) else max_thr_mad
    cap_applied <- is.finite(cap) && cap < max_thr_mad
  }

  # A zero MAD collapses the threshold onto the median, so ~half the group fails
  # on a metric with no spread. This is what a near-zero-inflated metric does:
  # percent.mt in single-NUCLEUS data has a median of ~0.03% and a MAD of 0, and
  # a 3-MAD ceiling then deletes any cell with any mitochondrial content at all.
  # The guardrails are what stop it -- for a high_bad metric, `hard_floor` means
  # "never be more aggressive than this". Flagged in `method` so a degenerate
  # threshold is visible in the CSV rather than only in the cell counts.
  list(
    min_thr     = as.numeric(min_thr),
    max_thr     = as.numeric(max_thr),
    max_thr_mad = as.numeric(max_thr_mad),
    median_raw  = .from_scaled(med_s, scale),
    mad_raw     = .from_scaled(med_s + mad_s, scale) - .from_scaled(med_s, scale),
    cap_applied = isTRUE(cap_applied),
    method      = if (mad_s == 0) "mad_zero_spread" else "mad"
  )
}


# ── The metric-configuration contract ────────────────────────────────────────

#' Validate a QC metric configuration
#'
#' A metric config is a named list; each element describes one `meta.data` column
#' with these fields: `direction` (`"low_bad"` computes a minimum, `"high_bad"` a
#' maximum; required), `scale` (`"log2"`, the default, applies MAD on
#' `log2(x + 1)`; `"raw"` for signed or percentage metrics), `n_mad` (MADs from
#' the median, default 3), and the optional guardrails `hard_floor` and
#' `hard_ceiling`, which a threshold may never cross.
#'
#' @param metrics Named list of metric configurations.
#' @param meta Optional cell metadata (data.frame or Seurat object). When
#'   supplied, every metric name must be a numeric column of it.
#' @param what Label used in error messages.
#' @return `TRUE`, invisibly. Errors with an informative message otherwise.
#' @examples
#' validate_qc_metrics(list(
#'   nCount_RNA = list(direction = "low_bad"),
#'   percent.mt = list(direction = "high_bad", scale = "raw", hard_floor = 5)
#' ))
#' @export
validate_qc_metrics <- function(metrics, meta = NULL, what = "metric config") {

  if (!is.list(metrics) || !length(metrics))
    stop(what, " must be a non-empty named list", call. = FALSE)

  nms <- names(metrics)
  if (is.null(nms) || any(!nzchar(trimws(nms))))
    stop(what, " must be fully named (one name per meta.data column)",
         call. = FALSE)
  if (anyDuplicated(nms))
    stop(what, " has duplicated metric name(s): ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "), call. = FALSE)

  for (m in nms) {
    cfg <- metrics[[m]]
    if (!is.list(cfg))
      stop(what, ": '", m, "' must be a list", call. = FALSE)

    dir <- cfg$direction
    if (!is.character(dir) || length(dir) != 1L ||
        !dir %in% c("low_bad", "high_bad"))
      stop(what, ": '", m, "' has direction = ",
           if (is.null(dir)) "NULL" else paste(dQuote(dir), collapse = ", "),
           "; must be exactly \"low_bad\" or \"high_bad\"", call. = FALSE)

    sc <- cfg$scale %||% "log2"
    if (!is.character(sc) || length(sc) != 1L || !sc %in% c("log2", "raw"))
      stop(what, ": '", m, "' has scale = ", paste(sc, collapse = ", "),
           "; must be \"log2\" or \"raw\"", call. = FALSE)

    nm <- cfg$n_mad %||% 3
    if (!is.numeric(nm) || length(nm) != 1L || !is.finite(nm) || nm <= 0)
      stop(what, ": '", m, "' has n_mad = ", paste(nm, collapse = ", "),
           "; must be one positive finite number", call. = FALSE)

    for (g in c("hard_floor", "hard_ceiling")) {
      v <- cfg[[g]]
      if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || !is.finite(v)))
        stop(what, ": '", m, "' has ", g, " = ", paste(v, collapse = ", "),
             "; must be NULL or one finite number", call. = FALSE)
    }
    if (!is.null(cfg$hard_floor) && !is.null(cfg$hard_ceiling) &&
        cfg$hard_floor > cfg$hard_ceiling)
      stop(what, ": '", m, "' has hard_floor (", cfg$hard_floor,
           ") above hard_ceiling (", cfg$hard_ceiling, ")", call. = FALSE)
  }

  if (!is.null(meta)) {
    meta <- .qc_meta(meta)
    miss <- setdiff(nms, colnames(meta))
    if (length(miss))
      stop(what, ": metric column(s) not in the cell metadata: ",
           paste(miss, collapse = ", "), call. = FALSE)
    nonnum <- nms[!vapply(nms, function(m) is.numeric(meta[[m]]), logical(1))]
    if (length(nonnum))
      stop(what, ": metric column(s) are not numeric: ",
           paste(nonnum, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}


# ── Tier 1 : flag_bad_qc_groups() ────────────────────────────────────────────

#' Tier 1: flag systematically poor-quality groups
#'
#' Each group (e.g. cluster) is summarised by its median for every metric; MAD is
#' then applied *across those group medians*. A group whose median is an outlier
#' relative to the other groups is flagged, which catches debris clusters that a
#' per-cell cutoff computed within the same cluster cannot see. A group that is a
#' high outlier for both an RNA and an ATAC count metric is tagged
#' `possible_doublet` rather than `low_qc`. With fewer than three groups nothing
#' is flagged on metric values.
#'
#' @param x Cell metadata (data.frame) or Seurat object.
#' @param group_col Metadata column defining the groups.
#' @param metrics Metric config; see [validate_qc_metrics()].
#' @param n_mad_across MADs from the median of the group medians.
#' @param min_group_size Groups smaller than this get `review_small_n`.
#' @param n_bad_for_removal Number of bad metrics required for `low_qc`; fewer
#'   (but at least one) gives `review`.
#' @param rna_count_pattern,atac_count_pattern Regexes identifying the count
#'   metrics used for the doublet tag.
#' @return A data.frame with columns `group`, `n_cells`, `status`, `flag_reason`
#'   and one `median_<metric>` column per metric.
#' @seealso [compute_qc_thresholds()], [apply_qc_thresholds()]
#' @examples
#' set.seed(1)
#' meta <- data.frame(cl = rep(1:5, each = 40),
#'                    nCount_RNA = rpois(200, 3000),
#'                    percent.mt = runif(200, 0, 5))
#' meta$nCount_RNA[meta$cl == 5] <- 200
#' metrics <- list(nCount_RNA = list(direction = "low_bad"),
#'                 percent.mt = list(direction = "high_bad", scale = "raw"))
#' flag_bad_qc_groups(meta, "cl", metrics)[, 1:4]
#' @export
flag_bad_qc_groups <- function(x,
                               group_col,
                               metrics,
                               n_mad_across       = 3,
                               min_group_size     = 20,
                               n_bad_for_removal  = 2,
                               rna_count_pattern  = "nCount_RNA|nFeature_RNA",
                               atac_count_pattern = "nCount_ATAC|nFeature_ATAC") {

  meta <- .qc_meta(x)
  validate_qc_metrics(metrics, meta, what = "flag_bad_qc_groups metrics")
  if (!group_col %in% colnames(meta))
    stop("group_col '", group_col, "' not found in the cell metadata",
         call. = FALSE)

  metric_names <- names(metrics)
  grp    <- as.character(meta[[group_col]])
  groups <- .sort_groups(grp)
  if (!length(groups))
    stop("group_col '", group_col, "' has no non-missing values", call. = FALSE)

  # Per-group medians, computed in each metric's own scale then returned raw.
  med_mat <- matrix(NA_real_, nrow = length(metric_names), ncol = length(groups),
                    dimnames = list(metric_names, groups))
  for (m in metric_names) {
    sc <- metrics[[m]]$scale %||% "log2"
    for (g in groups) {
      v <- .finite(meta[[m]][grp == g])
      if (length(v))
        med_mat[m, g] <- .from_scaled(stats::median(.to_scaled(v, sc)), sc)
    }
  }

  # Cross-group bounds. Fewer than three groups cannot support a MAD across
  # them, so nothing is flagged on metric values -- recorded, not silent.
  bounds <- lapply(metric_names, function(m) {
    sc   <- metrics[[m]]$scale %||% "log2"
    vals <- .finite(med_mat[m, ])
    if (length(vals) < 3) return(list(lo = -Inf, hi = Inf, n = length(vals)))
    vs <- .to_scaled(vals, sc)
    list(lo = .from_scaled(stats::median(vs) - n_mad_across * stats::mad(vs), sc),
         hi = .from_scaled(stats::median(vs) + n_mad_across * stats::mad(vs), sc),
         n  = length(vals))
  })
  names(bounds) <- metric_names

  n_by_group <- table(grp)

  rows <- lapply(groups, function(g) {

    n_cells  <- as.integer(n_by_group[[g]])
    bad_low  <- character(0)
    bad_high <- character(0)

    for (m in metric_names) {
      val <- med_mat[m, g]
      if (!is.finite(val)) next
      dir <- metrics[[m]]$direction
      if (dir == "low_bad"  && val < bounds[[m]]$lo) bad_low  <- c(bad_low,  m)
      if (dir == "high_bad" && val > bounds[[m]]$hi) bad_high <- c(bad_high, m)
    }

    all_bad <- c(bad_low, bad_high)
    dbl     <- any(grepl(rna_count_pattern,  bad_high)) &&
               any(grepl(atac_count_pattern, bad_high))

    status <-
      if (n_cells < min_group_size)               "review_small_n"
      else if (dbl)                               "possible_doublet"
      else if (length(all_bad) >= n_bad_for_removal) "low_qc"
      else if (length(all_bad) >= 1)              "review"
      else                                        "pass"

    out <- data.frame(group = g, n_cells = n_cells, status = status,
                      flag_reason = paste(all_bad, collapse = "; "),
                      stringsAsFactors = FALSE)
    for (m in metric_names) out[[paste0("median_", m)]] <- unname(med_mat[m, g])
    out
  })

  flags <- do.call(rbind, rows)
  stopifnot(all(flags$status %in% QC_GROUP_STATUSES))
  rownames(flags) <- NULL
  flags
}


# ── Tier 2 : compute_qc_thresholds() ─────────────────────────────────────────

#' Tier 2: per-group, per-metric MAD thresholds for individual cells
#'
#' Thresholds are computed only for `passing_groups`, so a cluster that tier 1
#' marked for removal never acquires bounds; [qc_fail_reasons()] then insists such
#' a cluster be removed rather than pass unfiltered.
#'
#' @param x Cell metadata or Seurat object.
#' @param split_by Metadata column defining the groups.
#' @param metrics Metric config; see [validate_qc_metrics()].
#' @param passing_groups Groups to compute for; `NULL` means all.
#' @param parent_col Optional coarser grouping. A group with fewer than
#'   `min_cells_per_group` cells borrows its parent's pooled distribution.
#' @param min_cells_per_group Cell floor for the parent fallback.
#' @param upper_quantile `NULL` (default) or a probability in (0, 1] capping every
#'   metric at that quantile of the group. This removes a fixed `1 - q` fraction
#'   of each group by construction.
#' @return A data.frame with columns `split_by`, `group`, `metric`, `min_thr`,
#'   `max_thr`, `max_thr_mad`, `median_raw`, `mad_raw`, `n_cells`, `n_used`,
#'   `cap_applied`, `upper_quantile` and `method`.
#' @seealso [flag_bad_qc_groups()], [qc_fail_reasons()]
#' @examples
#' set.seed(1)
#' meta <- data.frame(cl = rep(1:3, each = 50), nCount_RNA = rpois(150, 3000))
#' compute_qc_thresholds(meta, "cl", list(nCount_RNA = list(direction = "low_bad")))
#' @export
compute_qc_thresholds <- function(x,
                                  split_by,
                                  metrics,
                                  passing_groups      = NULL,
                                  parent_col          = NULL,
                                  min_cells_per_group = 20,
                                  upper_quantile      = NULL) {

  meta <- .qc_meta(x)
  validate_qc_metrics(metrics, meta, what = "compute_qc_thresholds metrics")
  if (!split_by %in% colnames(meta))
    stop("split_by '", split_by, "' not found in the cell metadata", call. = FALSE)
  if (!is.null(parent_col) && !parent_col %in% colnames(meta))
    stop("parent_col '", parent_col, "' not found in the cell metadata",
         call. = FALSE)
  if (!is.null(upper_quantile) &&
      (!is.numeric(upper_quantile) || length(upper_quantile) != 1L ||
       !is.finite(upper_quantile) || upper_quantile <= 0 || upper_quantile > 1))
    stop("upper_quantile must be NULL or one number in (0, 1]; got: ",
         paste(upper_quantile, collapse = ", "), call. = FALSE)

  grp    <- as.character(meta[[split_by]])
  groups <- .sort_groups(grp)
  if (!is.null(passing_groups))
    groups <- intersect(groups, as.character(passing_groups))

  rows <- list()
  for (g in groups) {
    idx     <- which(grp == g)
    n_cells <- length(idx)

    for (m in names(metrics)) {
      v      <- .finite(meta[[m]][idx])
      suffix <- ""

      if (length(v) < min_cells_per_group && !is.null(parent_col)) {
        parent  <- as.character(meta[[parent_col]])[idx][1]
        v_par   <- .finite(meta[[m]][as.character(meta[[parent_col]]) %in% parent])
        if (length(v_par) >= min_cells_per_group) {
          v <- v_par; suffix <- "_pooled_parent"
        } else {
          suffix <- "_too_few"
        }
      }

      r <- .compute_one(v, metrics[[m]], upper_quantile)
      rows[[length(rows) + 1L]] <- data.frame(
        split_by       = split_by,
        group          = g,
        metric         = m,
        min_thr        = r$min_thr,
        max_thr        = r$max_thr,
        max_thr_mad    = r$max_thr_mad,
        median_raw     = r$median_raw,
        mad_raw        = r$mad_raw,
        n_cells        = n_cells,
        n_used         = length(v),
        cap_applied    = r$cap_applied,
        upper_quantile = if (is.null(upper_quantile)) NA_real_ else upper_quantile,
        method         = paste0(r$method, suffix),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows))
    stop("compute_qc_thresholds produced no rows -- no group in '", split_by,
         "' survived `passing_groups`. Every cell would be unfiltered.",
         call. = FALSE)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# ── Tier 2 applied : qc_fail_reasons() ───────────────────────────────────────

.fail_reasons_one <- function(meta, thresholds, cluster_col, modality,
                              exempt_groups = character(0)) {

  req <- c("group", "metric", "min_thr", "max_thr")
  miss <- setdiff(req, colnames(thresholds))
  if (length(miss))
    stop("threshold table for '", modality, "' is missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  if (!cluster_col %in% colnames(meta))
    stop("cluster column '", cluster_col, "' not found in the cell metadata",
         call. = FALSE)

  metrics <- unique(thresholds$metric)
  miss_m  <- setdiff(metrics, colnames(meta))
  if (length(miss_m))
    stop("threshold table for '", modality, "' names metric(s) absent from the ",
         "cell metadata: ", paste(miss_m, collapse = ", "), call. = FALSE)

  grp <- as.character(meta[[cluster_col]])
  thresholds$group <- as.character(thresholds$group)

  # A cluster with no bounds and no removal flag would pass every cell
  # unfiltered. That is the defect this guard exists for.
  unbounded <- setdiff(unique(grp), c(thresholds$group, as.character(exempt_groups)))
  if (length(unbounded))
    stop("no QC thresholds for ", length(unbounded), " '", cluster_col,
         "' group(s), and they are not flagged for removal: ",
         paste(utils::head(unbounded, 8), collapse = ", "),
         "\n  Cells in these groups would pass QC unfiltered. Either compute ",
         "thresholds for them or flag them.", call. = FALSE)

  reasons <- rep("", nrow(meta))
  for (m in metrics) {
    sub   <- thresholds[thresholds$metric == m, , drop = FALSE]
    i     <- match(grp, sub$group)
    v     <- meta[[m]]
    vmin  <- sub$min_thr[i]
    vmax  <- sub$max_thr[i]

    fail_lo <- !is.na(vmin) & is.finite(vmin) & !is.na(v) & v < vmin
    fail_hi <- !is.na(vmax) & is.finite(vmax) & !is.na(v) & v > vmax

    reasons <- .paste_reason(reasons, ifelse(fail_lo, paste0(m, "<min"), ""))
    reasons <- .paste_reason(reasons, ifelse(fail_hi, paste0(m, ">max"), ""))
  }
  reasons
}

# Semicolon-join two reason vectors, treating "" as absent.
.paste_reason <- function(a, b) {
  ifelse(a == "", b, ifelse(b == "", a, paste0(a, ";", b)))
}

#' One QC verdict per cell, as a reason string
#'
#' Pure function: applies tier-1 group flags and tier-2 thresholds, and optionally
#' doublet calls. Errors if a group has neither thresholds nor a removal flag,
#' since its cells would otherwise pass unfiltered.
#'
#' @param x Cell metadata or Seurat object.
#' @param modality Named list, one entry per modality. Each entry is a list with
#'   `thresholds` (from [compute_qc_thresholds()]), `cluster_col`, and optionally
#'   `flagged_groups` (whole groups to remove). The name appears in the reason,
#'   e.g. `"low_qc_RNA_group"`.
#' @param remove_doublets Whether `doublet_col == "doublet"` fails a cell.
#' @param doublet_col Metadata column holding the scDblFinder class.
#' @return Character vector, one element per cell; `""` means the cell passes.
#' @examples
#' set.seed(1)
#' meta <- data.frame(cl = rep(1:3, each = 50), nCount_RNA = rpois(150, 3000))
#' thr <- compute_qc_thresholds(meta, "cl",
#'                              list(nCount_RNA = list(direction = "low_bad")))
#' table(qc_fail_reasons(meta, list(RNA = list(thresholds = thr, cluster_col = "cl")),
#'                       remove_doublets = FALSE))
#' @export
qc_fail_reasons <- function(x,
                            modality,
                            remove_doublets = TRUE,
                            doublet_col     = "scDblFinder.class") {

  meta <- .qc_meta(x)
  if (!is.list(modality) || !length(modality) || is.null(names(modality)))
    stop("`modality` must be a non-empty named list", call. = FALSE)

  reasons <- rep("", nrow(meta))

  for (nm in names(modality)) {
    spec <- modality[[nm]]
    for (f in c("thresholds", "cluster_col"))
      if (is.null(spec[[f]]))
        stop("modality '", nm, "' is missing `", f, "`", call. = FALSE)

    flagged <- as.character(spec$flagged_groups %||% character(0))
    if (length(flagged)) {
      in_bad  <- as.character(meta[[spec$cluster_col]]) %in% flagged
      reasons <- .paste_reason(reasons,
                               ifelse(in_bad, paste0("low_qc_", nm, "_group"), ""))
    }

    reasons <- .paste_reason(
      reasons,
      .fail_reasons_one(meta, spec$thresholds, spec$cluster_col, nm,
                        exempt_groups = flagged))
  }

  if (remove_doublets) {
    if (!doublet_col %in% colnames(meta))
      stop("remove_doublets = TRUE but '", doublet_col, "' is not in the cell ",
           "metadata. Run scDblFinder first, or pass remove_doublets = FALSE ",
           "deliberately.", call. = FALSE)
    is_dbl <- !is.na(meta[[doublet_col]]) & meta[[doublet_col]] == "doublet"
    reasons <- .paste_reason(reasons, ifelse(is_dbl, "doublet", ""))
  }

  reasons
}


#' Apply QC to a Seurat object: verdicts, subset and summaries
#'
#' Adds a `QC_fail_reason` metadata column, subsets to passing cells, and returns
#' summary tables. Refuses to return an empty object.
#'
#' @inheritParams qc_fail_reasons
#' @param seurat_obj Seurat object to filter.
#' @param return_removed_obj Also return the removed cells as an object.
#' @return A list with `obj` (filtered Seurat object), `fail_reason` (named
#'   character vector), `summary` (one-row data.frame of counts), and
#'   `removed_breakdown` (cells per reason); plus `removed_obj` if requested.
#' @seealso [qc_fail_reasons()]
#' @export
apply_qc_thresholds <- function(seurat_obj,
                                modality,
                                remove_doublets    = TRUE,
                                doublet_col        = "scDblFinder.class",
                                return_removed_obj = FALSE) {

  if (!isTRUE(methods::is(seurat_obj, "Seurat")))
    stop("apply_qc_thresholds() needs a Seurat object; use qc_fail_reasons() ",
         "for a plain metadata table", call. = FALSE)

  meta    <- seurat_obj@meta.data
  reasons <- qc_fail_reasons(meta, modality, remove_doublets, doublet_col)
  names(reasons) <- rownames(meta)

  seurat_obj$QC_fail_reason <- unname(reasons[colnames(seurat_obj)])

  keep <- rownames(meta)[reasons == ""]
  drop <- rownames(meta)[reasons != ""]
  if (!length(keep))
    stop("QC removed every cell -- refusing to return an empty object. ",
         "Inspect the threshold tables.", call. = FALSE)

  tab <- table(ifelse(reasons == "", "kept", reasons))
  breakdown <- data.frame(reason = names(tab), n_cells = as.integer(tab),
                          stringsAsFactors = FALSE)
  breakdown <- breakdown[order(-breakdown$n_cells), , drop = FALSE]
  rownames(breakdown) <- NULL

  out <- list(
    obj         = subset(seurat_obj, cells = keep),
    fail_reason = reasons,
    summary     = data.frame(
      total_cells   = nrow(meta),
      kept_cells    = length(keep),
      removed_cells = length(drop),
      pct_removed   = round(100 * length(drop) / nrow(meta), 2),
      stringsAsFactors = FALSE),
    removed_breakdown = breakdown
  )
  if (return_removed_obj && length(drop))
    out$removed_obj <- subset(seurat_obj, cells = drop)
  out
}
