# Synthetic multi-sample metadata with a tier-1/tier-2 verdict.
make_sample_meta <- function(seed, n_cl = 3, n_per = 60) {
  set.seed(seed)
  n <- n_cl * n_per
  m <- data.frame(cl = rep(seq_len(n_cl), each = n_per),
                  nCount_RNA   = rpois(n, 3000),
                  nFeature_RNA = rpois(n, 1500),
                  percent.mt   = runif(n, 0, 8),
                  row.names    = paste0("s", seed, "_c", seq_len(n)))
  thr <- compute_qc_thresholds(m, "cl", toy_metrics())
  m$QC_fail_reason <- qc_fail_reasons(
    m, list(RNA = list(thresholds = thr, cluster_col = "cl")),
    remove_doublets = FALSE)
  list(meta = m, thr = thr, flags = flag_bad_qc_groups(m, "cl", toy_metrics()))
}

export_toy <- function(qc_dir, ids = c("S1", "S2", "S10")) {
  for (i in seq_along(ids)) {
    s <- make_sample_meta(i)
    export_qc_sample(s$meta, ids[i], qc_dir, cluster_col = "cl",
                     thresholds = s$thr, group_flags = s$flags)
  }
}

test_that("export_qc_sample writes a small, complete file", {
  d <- withr::local_tempdir()
  s <- make_sample_meta(1)
  p <- suppressMessages(export_qc_sample(s$meta, "S1", d, cluster_col = "cl",
                                         thresholds = s$thr, group_flags = s$flags))
  expect_true(file.exists(p))
  obj <- readRDS(p)
  expect_equal(nrow(obj$cells), nrow(s$meta))
  expect_equal(obj$summary$n_cells, nrow(s$meta))
  expect_equal(obj$summary$n_fail, sum(s$meta$QC_fail_reason != ""))
  expect_equal(obj$summary$n_pass + obj$summary$n_fail, obj$summary$n_cells)
  expect_true(all(c("median_nCount_RNA", "median_nCount_RNA_pass") %in%
                    colnames(obj$summary)))
  expect_equal(unique(obj$thresholds$sample), "S1")
  expect_lt(file.size(p), 1e5)
})

test_that("export_qc_sample accepts a Seurat object and pre-QC input", {
  obj <- toy_seurat()
  obj$percent.mt <- 1
  d <- withr::local_tempdir()
  expect_warning(
    suppressMessages(export_qc_sample(obj, "T1", d, cluster_col = "cl",
                                      metrics = c("nCount_RNA", "percent.mt"))),
    "pre-QC")
  o <- readRDS(file.path(d, "cells", "T1.rds"))
  expect_true(all(is.na(o$cells$qc_pass)))
  expect_true(is.na(o$summary$n_fail))
})

test_that("export_qc_sample error paths", {
  d <- withr::local_tempdir()
  m <- make_sample_meta(1)$meta
  expect_error(export_qc_sample(m, "../x", d), "file name")
  expect_error(export_qc_sample(m, c("a", "b"), d), "file name")
  expect_error(export_qc_sample(m, "S1", d, metrics = "nope"), "not in the cell")
  expect_error(export_qc_sample(m, "S1", d, metrics = "cl2"), "not in the cell")
  m$chr <- "a"
  expect_error(export_qc_sample(m, "S1", d, metrics = "chr"), "not numeric")
  expect_error(export_qc_sample(m, "S1", d, cluster_col = "zzz"), "cluster_col")
  expect_error(export_qc_sample(m, "S1", d, thresholds = data.frame(a = 1)),
               "missing column")
  expect_error(export_qc_sample(m[0, ], "S1", d), "no cells")
  expect_error(export_qc_sample(1:3, "S1", d), "expected a data.frame")
})

test_that("build_qc_report writes combined CSVs and is incremental", {
  d <- withr::local_tempdir()
  suppressMessages(export_toy(d, c("S1", "S2")))
  res <- suppressMessages(build_qc_report(d, render = FALSE))
  expect_null(res$html)
  expect_equal(res$samples, c("S1", "S2"))
  expect_true(all(file.exists(res$csv)))
  expect_equal(nrow(read.csv(res$csv[["summary"]])), 2)
  bd <- read.csv(res$csv[["removal_breakdown"]])
  expect_equal(sum(bd$n_cells), sum(res$summary$n_cells))

  # add a sample -> rebuilt report includes it, natural sort order (S2 < S10)
  s <- make_sample_meta(3)
  suppressMessages(export_qc_sample(s$meta, "S10", d, cluster_col = "cl"))
  res2 <- suppressMessages(build_qc_report(d, render = FALSE))
  expect_equal(res2$samples, c("S1", "S2", "S10"))
  expect_equal(nrow(read.csv(res2$csv[["summary"]])), 3)

  res3 <- suppressMessages(build_qc_report(d, samples = "S2", render = FALSE))
  expect_equal(res3$samples, "S2")
})

test_that("build_qc_report error paths", {
  d <- withr::local_tempdir()
  expect_error(build_qc_report(file.path(d, "nope")), "does not exist")
  expect_error(build_qc_report(d), "no exported samples")
  suppressMessages(export_toy(d, "S1"))
  expect_error(build_qc_report(d, samples = "zz", render = FALSE), "not exported")
  expect_error(build_qc_report(d, max_cells_plot = 0, render = FALSE),
               "max_cells_plot")
  saveRDS(list(a = 1), file.path(d, "cells", "bad.rds"))
  expect_error(build_qc_report(d, render = FALSE), "not an export_qc_sample")
})

test_that("build_qc_report renders a self-contained HTML", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  d <- withr::local_tempdir()
  suppressMessages(export_toy(d))
  res <- suppressMessages(build_qc_report(d))
  expect_true(file.exists(res$html))
  html <- paste(readLines(res$html, warn = FALSE), collapse = "\n")
  expect_match(html, "Per-sample QC report")
  expect_match(html, "S10")
  expect_match(html, "data:image/png;base64")   # figures embedded
  # only the report dir is written to; qc_dir/cells untouched
  expect_setequal(list.files(d), c("cells", "report"))
})
