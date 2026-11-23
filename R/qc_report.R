# Lightweight per-sample QC export + combined HTML/CSV report.
#
# Design: each sample is exported once (export_qc_sample) to a small .rds that
# holds only a handful of per-cell QC columns -- never the Seurat object. The
# report (build_qc_report) is rebuilt from whatever .rds files are in the folder,
# so adding or re-thresholding a sample is: re-export that sample, rebuild.

.QC_EXPORT_FORMAT <- 1L

# Sample ids become file names, so keep them boring and path-safe.
.check_sample_id <- function(sample) {
  if (!is.character(sample) || length(sample) != 1L || is.na(sample) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", sample))
    stop("`sample` must be one string of letters, digits, '.', '_' or '-' ",
         "(it is used as a file name); got: ",
         paste(format(sample), collapse = ", "), call. = FALSE)
  invisible(sample)
}

# "S2" before "S10": order by the prefix, then by the trailing number.
.natural_order <- function(x) {
  num <- suppressWarnings(as.numeric(sub("^.*?([0-9]+)$", "\\1", x, perl = TRUE)))
  order(sub("[0-9]+$", "", x), num, x)
}

.check_table <- function(tab, required, what) {
  if (is.null(tab)) return(NULL)
  if (!is.data.frame(tab))
    stop("`", what, "` must be a data.frame (or NULL)", call. = FALSE)
  miss <- setdiff(required, colnames(tab))
  if (length(miss))
    stop("`", what, "` is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  tab
}

# One-row summary of a sample's cell table.
.qc_sample_summary <- function(cells, metrics) {
  has_qc <- any(!is.na(cells$qc_pass))
  n      <- nrow(cells)
  n_pass <- if (has_qc) sum(cells$qc_pass, na.rm = TRUE) else NA_integer_
  n_fail <- if (has_qc) sum(!cells$qc_pass, na.rm = TRUE) else NA_integer_
  has_dbl <- any(!is.na(cells$doublet_class))
  n_dbl  <- if (has_dbl) sum(cells$doublet_class == "doublet", na.rm = TRUE) else
    NA_integer_
  out <- data.frame(
    sample      = cells$sample[1],
    n_cells     = n,
    n_pass      = as.integer(n_pass),
    n_fail      = as.integer(n_fail),
    pct_fail    = round(100 * n_fail / n, 2),
    n_doublet   = as.integer(n_dbl),
    pct_doublet = round(100 * n_dbl / n, 2),
    stringsAsFactors = FALSE)
  med <- function(v) {
    v <- .finite(v)
    if (length(v)) stats::median(v) else NA_real_
  }
  for (m in metrics) out[[paste0("median_", m)]] <- med(cells[[m]])
  if (has_qc)
    for (m in metrics)
      out[[paste0("median_", m, "_pass")]] <- med(cells[[m]][cells$qc_pass %in% TRUE])
  out
}

#' Export a lightweight per-sample QC file
#'
#' Keeps only a few per-cell QC columns (counts, genes, percent mito, cluster,
#' doublet call, QC verdict) plus optional threshold and group-flag tables, and
#' writes them to one small `.rds` file. The Seurat object, assays and embeddings
#' are never stored, so the output is a few hundred KB per sample and cheap to
#' keep, download and re-read. Run it once per sample (re-running overwrites that
#' sample only), then call [build_qc_report()] to compare all samples.
#'
#' Pass the object or metadata **before** removing failing cells, with the
#' verdict in `fail_col` (e.g. from [qc_fail_reasons()]), so the report can show
#' what the thresholds remove. If `fail_col` is absent the sample is exported as
#' pre-QC (pass/fail columns are `NA`) with a warning.
#'
#' @param x Cell metadata (data.frame) or Seurat object for one sample.
#' @param sample Sample id; used as the file name, so letters, digits, `.`, `_`
#'   and `-` only.
#' @param qc_dir Folder holding all exported samples (created if needed). Files
#'   go to `qc_dir/cells/<sample>.rds`.
#' @param metrics Numeric metadata columns to keep and summarise.
#' @param cluster_col Optional cluster column, used for per-cluster plots.
#' @param fail_col Column with the per-cell QC verdict, `""` meaning pass (see
#'   [qc_fail_reasons()]).
#' @param doublet_col Optional column with the doublet class (`"doublet"` /
#'   `"singlet"`); ignored if absent.
#' @param thresholds Optional [compute_qc_thresholds()] output for this sample.
#' @param group_flags Optional [flag_bad_qc_groups()] output for this sample.
#' @return The path of the written `.rds` file, invisibly.
#' @seealso [build_qc_report()]
#' @examples
#' set.seed(1)
#' meta <- data.frame(cl = rep(1:3, each = 50),
#'                    nCount_RNA = rpois(150, 3000),
#'                    nFeature_RNA = rpois(150, 1500),
#'                    percent.mt = runif(150, 0, 8),
#'                    row.names = paste0("cell", 1:150))
#' metrics <- list(nCount_RNA = list(direction = "low_bad"),
#'                 percent.mt = list(direction = "high_bad", scale = "raw"))
#' thr <- compute_qc_thresholds(meta, "cl", metrics)
#' meta$QC_fail_reason <- qc_fail_reasons(
#'   meta, list(RNA = list(thresholds = thr, cluster_col = "cl")),
#'   remove_doublets = FALSE)
#' qc_dir <- file.path(tempdir(), "qc_example")
#' export_qc_sample(meta, "S1", qc_dir, cluster_col = "cl", thresholds = thr)
#' list.files(qc_dir, recursive = TRUE)
#' @export
export_qc_sample <- function(x,
                             sample,
                             qc_dir,
                             metrics     = c("nCount_RNA", "nFeature_RNA",
                                             "percent.mt"),
                             cluster_col = NULL,
                             fail_col    = "QC_fail_reason",
                             doublet_col = "scDblFinder.class",
                             thresholds  = NULL,
                             group_flags = NULL) {
  .check_sample_id(sample)
  if (!is.character(qc_dir) || length(qc_dir) != 1L || !nzchar(qc_dir))
    stop("`qc_dir` must be one folder path", call. = FALSE)
  meta <- .qc_meta(x)
  if (!nrow(meta)) stop("no cells in `x`", call. = FALSE)
  if (!is.character(metrics) || !length(metrics) || anyDuplicated(metrics))
    stop("`metrics` must be a non-empty character vector without duplicates",
         call. = FALSE)
  miss <- setdiff(metrics, colnames(meta))
  if (length(miss))
    stop("metric column(s) not in the cell metadata: ",
         paste(miss, collapse = ", "), call. = FALSE)
  nonnum <- metrics[!vapply(metrics, function(m) is.numeric(meta[[m]]), logical(1))]
  if (length(nonnum))
    stop("metric column(s) are not numeric: ", paste(nonnum, collapse = ", "),
         call. = FALSE)
  if (!is.null(cluster_col) &&
      (!is.character(cluster_col) || length(cluster_col) != 1L ||
       !cluster_col %in% colnames(meta)))
    stop("cluster_col not found in the cell metadata: ",
         paste(cluster_col, collapse = ", "), call. = FALSE)
  thresholds  <- .check_table(thresholds, c("group", "metric", "min_thr", "max_thr"),
                              "thresholds")
  group_flags <- .check_table(group_flags, c("group", "status"), "group_flags")

  reason <- if (fail_col %in% colnames(meta)) as.character(meta[[fail_col]]) else {
    warning("'", fail_col, "' not in the cell metadata; exporting '", sample,
            "' as pre-QC (no pass/fail information).", call. = FALSE)
    rep(NA_character_, nrow(meta))
  }
  cells <- data.frame(
    sample = sample,
    cell   = rownames(meta),
    cluster = if (is.null(cluster_col)) NA_character_ else
      as.character(meta[[cluster_col]]),
    stringsAsFactors = FALSE)
  for (m in metrics) cells[[m]] <- as.numeric(meta[[m]])
  cells$doublet_class  <- if (!is.null(doublet_col) && doublet_col %in% colnames(meta))
    as.character(meta[[doublet_col]]) else NA_character_
  cells$qc_fail_reason <- reason
  cells$qc_pass        <- ifelse(is.na(reason), NA, reason == "")

  if (!is.null(thresholds))  thresholds  <- cbind(sample = sample, thresholds)
  if (!is.null(group_flags)) group_flags <- cbind(sample = sample, group_flags)

  obj <- list(format      = .QC_EXPORT_FORMAT,
              sample      = sample,
              metrics     = metrics,
              cells       = cells,
              summary     = .qc_sample_summary(cells, metrics),
              thresholds  = thresholds,
              group_flags = group_flags,
              created     = Sys.time())

  dir.create(file.path(qc_dir, "cells"), recursive = TRUE, showWarnings = FALSE)
  path <- file.path(qc_dir, "cells", paste0(sample, ".rds"))
  saveRDS(obj, path)
  message("Saved: ", path, " (", format(structure(file.size(path),
                                                  class = "object_size"),
                                        units = "auto"), ")")
  invisible(path)
}

.read_qc_sample <- function(path) {
  obj <- readRDS(path)
  if (!is.list(obj) || !all(c("format", "sample", "metrics", "cells", "summary")
                            %in% names(obj)))
    stop("not an export_qc_sample() file: ", path, call. = FALSE)
  if (!identical(obj$format, .QC_EXPORT_FORMAT))
    stop("unsupported export format ", obj$format, " in ", path,
         "; re-export the sample", call. = FALSE)
  obj
}

#' Build a combined QC report (HTML + CSV) from exported samples
#'
#' Reads every `qc_dir/cells/*.rds` written by [export_qc_sample()], writes
#' combined CSV tables and renders one self-contained HTML report comparing the
#' samples (summary table, removal reasons, before/after violins, count-vs-gene
#' scatter, thresholds and per-sample cluster plots). Only the small exported
#' files are read, so it is cheap to rerun whenever a sample is added or
#' re-thresholded. The CSVs are written before rendering, so they exist even if
#' rendering fails.
#'
#' Rendering needs the `rmarkdown` and `knitr` packages and pandoc. Use
#' `render = FALSE` to write only the CSVs.
#'
#' @param qc_dir Folder passed to [export_qc_sample()].
#' @param output_dir Where to write the report and CSVs; defaults to
#'   `qc_dir/report`.
#' @param samples Optional character vector restricting the report to these
#'   samples; default is every exported sample.
#' @param report_file Name of the HTML file.
#' @param title Report title.
#' @param render Render the HTML report (`TRUE`) or write only the CSVs.
#' @param max_cells_plot Cells per sample drawn in the scatter plot (evenly
#'   thinned, no random numbers used), to keep the HTML small.
#' @return A list with `html` (path or `NULL`), `csv` (named character vector of
#'   paths), `samples` and `summary` (the combined summary data.frame), invisibly.
#' @seealso [export_qc_sample()]
#' @examples
#' \donttest{
#' set.seed(1)
#' qc_dir <- file.path(tempdir(), "qc_report_example")
#' for (s in c("S1", "S2")) {
#'   meta <- data.frame(nCount_RNA = rpois(200, 3000),
#'                      nFeature_RNA = rpois(200, 1500),
#'                      percent.mt = runif(200, 0, 8),
#'                      QC_fail_reason = "",
#'                      row.names = paste0(s, "_cell", 1:200))
#'   meta$QC_fail_reason[meta$percent.mt > 6] <- "percent.mt>max"
#'   export_qc_sample(meta, s, qc_dir)
#' }
#' if (requireNamespace("rmarkdown", quietly = TRUE) &&
#'     rmarkdown::pandoc_available()) {
#'   res <- build_qc_report(qc_dir)
#'   res$summary[, 1:5]
#' }
#' }
#' @export
build_qc_report <- function(qc_dir,
                            output_dir     = file.path(qc_dir, "report"),
                            samples        = NULL,
                            report_file    = "qc_report.html",
                            title          = "Per-sample QC report",
                            render         = TRUE,
                            max_cells_plot = 5000) {
  if (!is.character(qc_dir) || length(qc_dir) != 1L || !dir.exists(qc_dir))
    stop("`qc_dir` does not exist: ", paste(qc_dir, collapse = ", "), call. = FALSE)
  if (!is.numeric(max_cells_plot) || length(max_cells_plot) != 1L ||
      !is.finite(max_cells_plot) || max_cells_plot < 1)
    stop("`max_cells_plot` must be one positive number", call. = FALSE)
  files <- list.files(file.path(qc_dir, "cells"), pattern = "\\.rds$",
                      full.names = TRUE)
  if (!length(files))
    stop("no exported samples in ", file.path(qc_dir, "cells"),
         "; run export_qc_sample() first", call. = FALSE)
  names(files) <- sub("\\.rds$", "", basename(files))
  if (!is.null(samples)) {
    miss <- setdiff(samples, names(files))
    if (length(miss))
      stop("sample(s) not exported: ", paste(miss, collapse = ", "),
           "\n  available: ", paste(names(files), collapse = ", "), call. = FALSE)
    files <- files[samples]
  }
  files <- files[.natural_order(names(files))]
  objs  <- lapply(files, .read_qc_sample)
  bad <- names(files)[vapply(names(files), function(n) objs[[n]]$sample != n,
                             logical(1))]
  if (length(bad))
    stop("file name and stored sample id disagree for: ",
         paste(bad, collapse = ", "), call. = FALSE)

  bind <- function(field) {
    parts <- Filter(Negate(is.null), lapply(objs, `[[`, field))
    if (length(parts)) {
      out <- dplyr::bind_rows(parts)
      rownames(out) <- NULL
      as.data.frame(out)
    } else NULL
  }
  cells   <- bind("cells")
  summ    <- bind("summary")
  thr     <- bind("thresholds")
  flags   <- bind("group_flags")
  metrics <- unique(unlist(lapply(objs, `[[`, "metrics")))

  # Removal reasons per sample (only where pass/fail is known).
  known <- cells[!is.na(cells$qc_pass), , drop = FALSE]
  breakdown <- NULL
  if (nrow(known)) {
    r <- ifelse(known$qc_pass, "pass", known$qc_fail_reason)
    tab <- as.data.frame(table(sample = known$sample, reason = r),
                         responseName = "n_cells", stringsAsFactors = FALSE)
    tab <- tab[tab$n_cells > 0, , drop = FALSE]
    tab$pct_of_sample <- round(100 * tab$n_cells /
                                 as.numeric(table(known$sample)[tab$sample]), 2)
    breakdown <- tab[order(tab$sample, -tab$n_cells), , drop = FALSE]
    rownames(breakdown) <- NULL
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  csv <- c(summary = file.path(output_dir, "qc_summary_all_samples.csv"))
  utils::write.csv(summ, csv[["summary"]], row.names = FALSE)
  if (!is.null(breakdown)) {
    csv[["removal_breakdown"]] <- file.path(output_dir, "qc_removal_breakdown.csv")
    utils::write.csv(breakdown, csv[["removal_breakdown"]], row.names = FALSE)
  }
  if (!is.null(thr)) {
    csv[["thresholds"]] <- file.path(output_dir, "qc_thresholds_all_samples.csv")
    utils::write.csv(thr, csv[["thresholds"]], row.names = FALSE)
  }
  if (!is.null(flags)) {
    csv[["group_flags"]] <- file.path(output_dir, "qc_group_flags_all_samples.csv")
    utils::write.csv(flags, csv[["group_flags"]], row.names = FALSE)
  }
  message("Wrote ", length(csv), " CSV file(s) to ", output_dir)

  html <- NULL
  if (render) {
    for (pkg in c("rmarkdown", "knitr"))
      if (!requireNamespace(pkg, quietly = TRUE))
        stop("package '", pkg, "' is needed to render the report ",
             "(CSV files were already written). Install it or use render = FALSE.",
             call. = FALSE)
    if (!rmarkdown::pandoc_available())
      stop("pandoc is needed to render the HTML report (CSV files were already ",
           "written). Install pandoc, or use render = FALSE.", call. = FALSE)
    tpl <- system.file("rmarkdown", "qc_report.Rmd", package = "seurassist",
                       mustWork = TRUE)
    work <- tempfile("qc_report_")
    dir.create(work)
    on.exit(unlink(work, recursive = TRUE), add = TRUE)
    file.copy(tpl, work)
    qc <- list(cells = cells, summary = summ, thresholds = thr, flags = flags,
               breakdown = breakdown, metrics = metrics,
               samples = names(files), max_cells_plot = max_cells_plot)
    html <- rmarkdown::render(
      file.path(work, "qc_report.Rmd"),
      output_file     = report_file,
      output_dir      = normalizePath(output_dir),
      intermediates_dir = work,
      knit_root_dir   = work,
      params          = list(qc = qc, title = title),
      envir           = new.env(parent = globalenv()),
      quiet           = TRUE)
    message("Wrote report: ", html)
  }
  invisible(list(html = html, csv = csv, samples = names(files), summary = summ))
}
