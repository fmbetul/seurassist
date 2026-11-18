test_that("qc_violin_panel and qc_violin_stack return ggplots", {
  m <- toy_meta()
  expect_s3_class(qc_violin_panel(m, "cl", "nCount_RNA"), "ggplot")
  thr <- compute_qc_thresholds(m, "cl", toy_metrics())
  fl  <- flag_bad_qc_groups(m, "cl", toy_metrics())
  expect_s3_class(qc_violin_panel(m, "cl", "nCount_RNA", thr, fl,
                                  show_legend = TRUE), "ggplot")
  p <- qc_violin_stack(m, "cl", c("nCount_RNA", "percent.mt"))
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork") ||
              inherits(p, "gg"))
})

test_that("qc_status_umap and labels work", {
  set.seed(1)
  emb <- matrix(rnorm(400), 200, 2, dimnames = list(NULL, c("UMAP_1", "UMAP_2")))
  st <- factor(sample(c("pass", "fail_RNA"), 200, TRUE, c(.9, .1)),
               levels = c("pass", "fail_RNA", "fail_ATAC", "fail_both", "doublet"))
  expect_s3_class(qc_status_umap(emb, st, title = "t"), "ggplot")
  expect_equal(qc_metric_label("nCount_RNA"), "RNA UMIs")
  expect_s3_class(theme_qc(), "theme")
})
