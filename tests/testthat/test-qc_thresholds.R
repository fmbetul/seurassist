test_that("validate_qc_metrics accepts a good config and rejects bad ones", {
  expect_true(validate_qc_metrics(toy_metrics(), toy_meta()))
  expect_error(validate_qc_metrics(list()), "non-empty")
  expect_error(validate_qc_metrics(list(list(direction = "low_bad"))), "named")
  expect_error(validate_qc_metrics(list(a = list(direction = "up"))), "direction")
  expect_error(validate_qc_metrics(list(a = list(direction = "low_bad",
                                                 scale = "ln"))), "scale")
  expect_error(validate_qc_metrics(list(a = list(direction = "low_bad",
                                                 n_mad = -1))), "n_mad")
  expect_error(validate_qc_metrics(list(a = list(direction = "low_bad",
                                                 hard_floor = 5,
                                                 hard_ceiling = 1))),
               "hard_floor")
  expect_error(validate_qc_metrics(list(nope = list(direction = "low_bad")),
                                   toy_meta()), "not in the cell metadata")
})

test_that("validate_qc_remove_statuses catches typos", {
  expect_identical(validate_qc_remove_statuses(c("low_qc", "low_qc")), "low_qc")
  expect_error(validate_qc_remove_statuses("low_quality"), "low_quality")
  expect_error(validate_qc_remove_statuses(1), "character")
})

test_that("flag_bad_qc_groups flags the low-count cluster only", {
  fl <- flag_bad_qc_groups(toy_meta(), "cl", toy_metrics())
  expect_equal(nrow(fl), 5)
  expect_true(all(c("group", "n_cells", "status", "flag_reason") %in% names(fl)))
  expect_equal(fl$status[fl$group == "5"] != "pass", TRUE)
  expect_true(all(fl$status[fl$group != "5"] == "pass"))
})

test_that("flag_bad_qc_groups handles small groups and bad input", {
  m <- toy_meta()
  fl <- flag_bad_qc_groups(m, "cl", toy_metrics(), min_group_size = 100)
  expect_true(all(fl$status == "review_small_n"))
  expect_error(flag_bad_qc_groups(m, "missing", toy_metrics()), "not found")
})

test_that("compute_qc_thresholds returns one row per group x metric", {
  thr <- compute_qc_thresholds(toy_meta(), "cl", toy_metrics())
  expect_equal(nrow(thr), 5 * 2)
  expect_true(all(c("group", "metric", "min_thr", "max_thr") %in% names(thr)))
  expect_true(all(is.finite(thr$min_thr[thr$metric == "nCount_RNA"])))

  thr2 <- compute_qc_thresholds(toy_meta(), "cl", toy_metrics(),
                                passing_groups = c(1, 2))
  expect_setequal(unique(thr2$group), c("1", "2"))
})

test_that("compute_qc_thresholds errors on bad input", {
  m <- toy_meta()
  expect_error(compute_qc_thresholds(m, "nope", toy_metrics()), "not found")
  expect_error(compute_qc_thresholds(m, "cl", toy_metrics(),
                                     upper_quantile = 2), "upper_quantile")
  expect_error(compute_qc_thresholds(m, "cl", toy_metrics(),
                                     passing_groups = "zzz"), "no rows")
})

test_that("qc_fail_reasons marks flagged groups and low cells", {
  m <- toy_meta()
  thr <- compute_qc_thresholds(m, "cl", toy_metrics(), passing_groups = 1:4)
  r <- qc_fail_reasons(m, list(RNA = list(thresholds = thr, cluster_col = "cl",
                                          flagged_groups = "5")),
                       remove_doublets = FALSE)
  expect_length(r, nrow(m))
  expect_true(all(grepl("low_qc_RNA_group", r[m$cl == 5])))
  expect_true(mean(r[m$cl != 5] == "") > 0.9)
})

test_that("qc_fail_reasons validates modality and doublet column", {
  m <- toy_meta()
  expect_error(qc_fail_reasons(m, list()), "non-empty")
  expect_error(qc_fail_reasons(m, list(RNA = list(cluster_col = "cl")),
                               remove_doublets = FALSE), "thresholds")
  thr <- compute_qc_thresholds(m, "cl", toy_metrics())
  expect_error(qc_fail_reasons(m, list(RNA = list(thresholds = thr,
                                                  cluster_col = "cl"))),
               "scDblFinder")
  m$scDblFinder.class <- rep(c("singlet", "doublet"), length.out = nrow(m))
  r <- qc_fail_reasons(m, list(RNA = list(thresholds = thr, cluster_col = "cl")))
  expect_true(all(grepl("doublet", r[m$scDblFinder.class == "doublet"])))
})

test_that("apply_qc_thresholds subsets a Seurat object", {
  obj <- toy_seurat()
  obj$nCount_RNA[1:5] <- 1
  thr <- compute_qc_thresholds(obj, "cl",
                               list(nCount_RNA = list(direction = "low_bad")))
  res <- apply_qc_thresholds(obj,
                             list(RNA = list(thresholds = thr, cluster_col = "cl")),
                             remove_doublets = FALSE,
                             return_removed_obj = TRUE)
  expect_equal(res$summary$total_cells, 60)
  expect_equal(res$summary$kept_cells + res$summary$removed_cells, 60)
  expect_equal(ncol(res$obj), res$summary$kept_cells)
  expect_length(res$fail_reason, 60)
  expect_true(all(c("reason", "n_cells") %in% names(res$removed_breakdown)))
  expect_equal(ncol(res$removed_obj), res$summary$removed_cells)
  expect_error(apply_qc_thresholds(toy_meta(), list()), "Seurat object")
})
