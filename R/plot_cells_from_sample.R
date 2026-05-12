#' Plot cells of a specified sample and cluster type
#'
#' This function highlights cells from a specified cluster in a given sample
#' on the UMAP of an integrated Seurat object.
#'
#' @param sample A Seurat object that contains the sample data.
#' @param ClusterName A string specifying the cell type or cluster name to highlight.
#' @param integrated An integrated Seurat object containing merged data from multiple samples.
#' @param cluster_col A string specifying the metadata column that holds cluster or cell type
#'   labels. Default is "celltype".
#'
#' @return The UMAP \code{ggplot} object, invisibly. Also prints the plot to the console.
#'
#' @seealso \code{\link[Seurat]{DimPlot}}
#'
#' @importFrom Seurat DimPlot
#'
#' @export
plot_cells_from_sample <- function(sample, ClusterName, integrated, cluster_col = "celltype") {

  sample_name <- as.character(substitute(sample))

  if (!cluster_col %in% colnames(sample[[]])) {
    stop("Column '", cluster_col, "' not found in sample metadata.")
  }

  # Get cell barcodes for the specified cluster
  meta <- sample[[]][, cluster_col, drop = FALSE]
  cell.barcodes <- rownames(meta[meta[[cluster_col]] == ClusterName, , drop = FALSE])

  if (length(cell.barcodes) == 0) {
    warning("No cells found for cluster '", ClusterName, "' in sample.")
  }

  # Add sample prefix to match barcode naming in the integrated object
  cell.barcodes <- paste0(sample_name, "_", cell.barcodes)

  p <- DimPlot(integrated, group.by = cluster_col, reduction = "umap",
               label = TRUE, cells.highlight = cell.barcodes)
  print(p)
  invisible(p)
}
