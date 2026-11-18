#' Reorder a cluster column so numeric-looking levels sort numerically
#'
#' Lexical factor order puts cluster 10 between 1 and 2. This reorders a factor
#' metadata column numerically. Non-factor columns, and factors with any
#' non-numeric level, are returned untouched, so it is safe to apply across a
#' whole metadata table.
#'
#' @param obj Seurat object.
#' @param col Metadata column name.
#' @return The Seurat object with `col` re-levelled.
#' @seealso [sort_all_clusters()], [sorted_levels()]
#' @examples
#' \dontrun{
#' obj <- sort_clusters(obj, "seurat_clusters")
#' }
#' @export
sort_clusters <- function(obj, col) {
  if (!col %in% colnames(obj@meta.data))
    stop("column '", col, "' not found in the object metadata", call. = FALSE)
  lvls <- levels(obj@meta.data[[col]])
  if (is.null(lvls) || anyNA(suppressWarnings(as.numeric(lvls)))) return(obj)
  obj@meta.data[[col]] <- factor(obj@meta.data[[col]],
                                 levels = lvls[order(as.numeric(lvls))])
  obj
}

#' Apply [sort_clusters()] to every metadata column matching a pattern
#'
#' @param obj Seurat object.
#' @param pattern Regular expression matched against metadata column names.
#' @return The Seurat object with all matching columns re-levelled.
#' @examples
#' \dontrun{
#' obj <- sort_all_clusters(obj, pattern = "_clusters_res")
#' }
#' @export
sort_all_clusters <- function(obj, pattern = "_clusters_res") {
  for (cl in grep(pattern, colnames(obj@meta.data), value = TRUE))
    obj <- sort_clusters(obj, cl)
  obj
}
