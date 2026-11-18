test_that("count_percent_by counts, sorts and errors", {
  obj <- toy_seurat()
  obj$ct <- rep(c("A", "B", "C"), times = c(30, 20, 10))
  df <- count_percent_by(obj, "ct")
  expect_equal(names(df), c("CellType", "Count", "Percent"))
  expect_equal(as.character(df$CellType), c("A", "B", "C"))
  expect_equal(df$Count, c(30, 20, 10))
  expect_equal(df$Percent, c(50, 33.3, 16.7))
  expect_equal(names(count_percent_by(obj, "ct", colname = "X"))[1], "X")
  expect_error(count_percent_by(obj, "nope"), "not found")
})

test_that("sorted_levels is numeric-aware", {
  expect_equal(sorted_levels(c("10", "2", "1")), c("1", "2", "10"))
  expect_equal(sorted_levels(c("b", "a")), c("a", "b"))
  expect_equal(sorted_levels(factor(c("b", "a"), levels = c("b", "a"))),
               c("b", "a"))
})

test_that("sort_clusters re-levels numeric factors only", {
  obj <- toy_seurat()
  obj$f <- factor(obj$cl)                       # lexical: 1, 10, 2
  expect_equal(levels(sort_clusters(obj, "f")$f), c("1", "2", "10"))
  obj$g <- factor(rep(c("x", "y", "z"), each = 20))
  expect_equal(levels(sort_clusters(obj, "g")$g), c("x", "y", "z"))
  expect_error(sort_clusters(obj, "nope"), "not found")
})

test_that("sort_all_clusters applies to matching columns", {
  obj <- toy_seurat()
  obj$a_clusters_res0.5 <- factor(obj$cl)
  obj$other <- factor(obj$cl)
  out <- sort_all_clusters(obj)
  expect_equal(levels(out$a_clusters_res0.5), c("1", "2", "10"))
  expect_equal(levels(out$other), levels(obj$other))
})

test_that("choose_dims returns a sensible range", {
  sd <- c(10, 8, 5, 3, 2, 1, rep(0.5, 20))
  d <- choose_dims(sd)
  expect_true(is.numeric(d) && d[1] == 1 && length(d) >= 10)
  d2 <- choose_dims(c(sd, rep(0.5, 20)), start = 2, cap_dims = 12)
  expect_equal(d2[1], 2)
  expect_lte(max(d2), 12)
  expect_error(choose_dims("a"), "numeric")
  expect_error(choose_dims(1), "at least 2")
  expect_error(choose_dims(c(1, -1)), "non-negative")
})
