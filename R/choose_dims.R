# choose_dims(): data-driven number of PCA/LSI dimensions.
# Adapted from the 28-novo-in-vivo pipeline.

#' Choose the number of PCA/LSI dimensions from component standard deviations
#'
#' A data-driven replacement for a fixed `1:30`. Takes the smaller of two
#' criteria: `co1`, the first dimension where cumulative variance passes
#' `var_threshold` and the dimension itself contributes less than `tail_contrib`
#' percent; and `co2`, the last dimension before the per-dimension change in
#' variance drops below `min_contrib`. For LSI use `start = 2`, since component 1
#' tracks sequencing depth.
#'
#' @param stdev Numeric vector of component standard deviations, e.g.
#'   `Seurat::Stdev(obj, "pca")`.
#' @param var_threshold Cumulative percent variance target.
#' @param min_contrib Per-dimension percent change floor defining the elbow.
#' @param tail_contrib Per-dimension percent ceiling for the cumulative criterion.
#' @param start First component considered; 2 for LSI.
#' @param floor_dims Never select fewer than this many dimensions.
#' @param cap_dims Never select more than this many (`NULL` = no cap).
#' @return Integer vector of dimensions, e.g. `1:k` for PCA or `2:k` for LSI.
#' @examples
#' choose_dims(c(10, 8, 5, 3, 2, 1, rep(0.5, 20)))
#' choose_dims(c(10, 8, 5, 3, 2, 1, rep(0.5, 40)), start = 2, cap_dims = 12)
#' @export
choose_dims <- function(stdev,
                        var_threshold = 90,
                        min_contrib   = 0.1,
                        tail_contrib  = 5,
                        start         = 1,
                        floor_dims    = 10,
                        cap_dims      = NULL) {

  if (!is.numeric(stdev) || length(stdev) < 2)
    stop("`stdev` must be a numeric vector of at least 2 components",
         call. = FALSE)
  if (any(!is.finite(stdev)) || any(stdev < 0))
    stop("`stdev` must be finite and non-negative", call. = FALSE)
  start <- as.integer(start)
  if (start < 1 || start > length(stdev) - 1)
    stop("`start` must be between 1 and length(stdev) - 1", call. = FALSE)

  s    <- stdev[start:length(stdev)]
  v    <- s^2
  if (sum(v) <= 0) stop("all components have zero variance", call. = FALSE)
  pct  <- v / sum(v) * 100
  cumu <- cumsum(pct)

  co1 <- which(cumu > var_threshold & pct < tail_contrib)[1]

  dpc <- pct[-length(pct)] - pct[-1]
  co2 <- sort(which(dpc > min_contrib), decreasing = TRUE)[1] + 1

  k <- suppressWarnings(min(c(co1, co2), na.rm = TRUE))
  if (!is.finite(k)) k <- length(s)

  cap <- if (is.null(cap_dims)) length(s) else cap_dims
  k   <- min(max(k, floor_dims), cap, length(s))

  start:(start + k - 1L)
}
