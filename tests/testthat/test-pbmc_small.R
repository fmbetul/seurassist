# Functions that write files or need a full Seurat object, exercised on the
# 80-cell pbmc_small object that ships with SeuratObject. Each test works in a
# temporary directory so no files leak into the package.

test_that("make_DotPlot returns a plot and per-panel gene counts", {
  obj <- pbmc_toy()
  g <- pbmc_genes(4)
  res <- suppressWarnings(suppressMessages(make_DotPlot(
    obj,
    list(A = g[1:2], B = c(g[2:3], "NOT_A_GENE"), C = "NOT_A_GENE"),
    group_by = "groups")))   # warns: only 2 groups to scale over
  expect_s3_class(res$plot, "ggplot")
  # g[2] is de-duplicated out of B; the absent gene and empty panel C are dropped
  expect_equal(res$groups_kept, c(A = 2L, B = 1L))
  expect_equal(SeuratObject::DefaultAssay(obj), "RNA")
})

test_that("make_DotPlot validates its input", {
  obj <- pbmc_toy()
  g <- pbmc_genes(2)
  expect_error(make_DotPlot(obj, list(A = g), assay = "nope"), "not found")
  expect_error(make_DotPlot(obj, list(A = g), group_by = "nope"),
               "Metadata column")
  expect_error(suppressMessages(make_DotPlot(obj, list(A = "NOT_A_GENE"))),
               "No genes remain")
  expect_error(make_DotPlot(obj, list()))
})

test_that("generate_scaled_FPs writes separate and combined PDFs", {
  withr::local_dir(withr::local_tempdir())
  objs <- list(a = pbmc_toy(), b = pbmc_toy())
  g <- pbmc_genes(2)

  quietly <- function(x) invisible(utils::capture.output(suppressMessages(x)))
  quietly(generate_scaled_FPs(objs, g, filename = "t", separate = TRUE))
  expect_true(all(file.exists(c("t_a_FPs.pdf", "t_b_FPs.pdf"))))

  quietly(generate_scaled_FPs(objs, g, filename = "t", separate = FALSE))
  expect_true(file.exists("t_FPs.pdf"))

  expect_error(generate_scaled_FPs(unname(objs), g, "t"), "named list")
  # FeaturePlot also warns about the missing gene; only check ours is raised
  msgs <- character()
  withCallingHandlers(
    quietly(generate_scaled_FPs(objs, c(g[1], "NOT_A_GENE"), filename = "w")),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("not found in any Seurat object", msgs)))
})

test_that("generate_AllFPs and generate_allVlns write PDFs", {
  withr::local_dir(withr::local_tempdir())
  pbmc <- pbmc_toy()
  gl <- list(first = pbmc_genes(3), none = "NOT_A_GENE")

  out <- utils::capture.output(
    expect_message(generate_AllFPs(pbmc, gl), "Skipping 'none'"))
  expect_true(file.exists("output/pbmc_AllFPs.pdf"))

  out <- utils::capture.output(
    expect_message(generate_allVlns(pbmc, gl), "Skipping 'none'"))
  expect_true(file.exists("output/pbmc_allVlns.pdf"))
})

test_that("save_qc_file writes a PDF and optional marker CSVs", {
  withr::local_dir(withr::local_tempdir())
  pbmc <- pbmc_toy()

  suppressMessages(suppressWarnings(save_qc_file(pbmc, version = "v0")))
  expect_true(file.exists("output/pbmc_qc_file_v0.pdf"))

  expect_warning(
    suppressMessages(save_qc_file(pbmc, version = "v1", qc.plot = FALSE,
                                  csv = TRUE)),
    "csv = TRUE has no effect")

  suppressMessages(suppressWarnings(
    save_qc_file(pbmc, version = "v2", qc.plot = FALSE, Heatmap = TRUE,
                 csv = TRUE)))
  expect_true(file.exists("output/pbmc_MarkersPerCluster_v2.csv"))
  expect_true(file.exists("output/pbmc_MarkersPerCluster_top10_v2.csv"))
})

test_that("plot_cells_from_sample and plot_cells_from_integrated highlight cells", {
  withr::local_dir(withr::local_tempdir())
  s1 <- pbmc_toy()
  integrated <- SeuratObject::RenameCells(s1, add.cell.id = "s1")

  p1 <- plot_cells_from_sample(s1, "g1", integrated, cluster_col = "groups")
  expect_s3_class(p1, "ggplot")

  p2 <- plot_cells_from_integrated(integrated, "g1", s1, cluster_col = "groups")
  expect_s3_class(p2, "ggplot")

  expect_error(plot_cells_from_sample(s1, "g1", integrated, "nope"),
               "not found in sample")
  expect_error(plot_cells_from_integrated(integrated, "g1", s1, "nope"),
               "not found in integrated")
  expect_warning(
    plot_cells_from_sample(s1, "no_such_cluster", integrated,
                           cluster_col = "groups"),
    "No cells found")
})
