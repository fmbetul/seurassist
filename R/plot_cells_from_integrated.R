#' Plot specific cells from an integrated Seurat object
#'
#' This function highlights cells from a specified cluster and sample on the
#' UMAP of a per-sample Seurat object, using barcode information from an
#' integrated object.
#'
#' @param integrated A Seurat object. An integrated or merged object that includes
#'   data from multiple samples.
#' @param ClusterName A character string specifying the name of the cell type or
#'   cluster whose cells are to be visualised.
#' @param sample A Seurat object containing scRNA-seq data for a specific sample.
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
plot_cells_from_integrated <- function(integrated, ClusterName, sample, cluster_col = "celltype") {

  sample_name <- as.character(substitute(sample))

  if (!cluster_col %in% colnames(integrated[[]])) {
    stop("Column '", cluster_col, "' not found in integrated metadata.")
  }

  # Get barcodes for the target cluster in the integrated object
  meta <- integrated[[]][, cluster_col, drop = FALSE]
  cell.barcodes <- rownames(meta[meta[[cluster_col]] == ClusterName, , drop = FALSE])

  # Retain only barcodes belonging to the target sample
  cell.barcodes <- grep(paste0("^", sample_name, "_"), cell.barcodes, value = TRUE)

  # Strip the sample prefix for lookup in the per-sample object
  cell.barcodes <- sub(paste0("^", sample_name, "_"), "", cell.barcodes)

  if (length(cell.barcodes) == 0) {
    warning("No cells found for cluster '", ClusterName, "' in sample '", sample_name, "'.")
  }

  DefaultAssay(sample) <- "RNA"

  if (cluster_col %in% colnames(sample[[]])) {
    p <- DimPlot(sample, group.by = cluster_col, reduction = "umap",
                 label = TRUE, cells.highlight = cell.barcodes)
  } else {
    p <- DimPlot(sample, label = TRUE, cells.highlight = cell.barcodes)
  }

  print(p)
  invisible(p)
}
