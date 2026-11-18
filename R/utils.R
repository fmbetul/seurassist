# Internal helpers shared across seurassist.

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Numeric-aware ordering of cluster labels
#'
#' An existing factor keeps its level order; character/numeric vectors sort
#' numerically when every label looks numeric ("2" before "10"), else lexically.
#'
#' @param x Character, numeric or factor vector.
#' @return Character vector of ordered unique labels.
#' @export
sorted_levels <- function(x) {
  if (is.factor(x)) return(intersect(levels(x), unique(as.character(x))))
  u <- unique(as.character(x))
  u <- u[!is.na(u)]
  n <- suppressWarnings(as.numeric(u))
  if (anyNA(n)) sort(u) else u[order(n)]
}
