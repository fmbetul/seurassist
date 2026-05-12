#' Generate Scaled Feature Plots
#'
#' This function generates feature plots for a list of genes in multiple Seurat objects,
#' scaling the colour axis by the maximum expression value of each gene across all objects.
#' The generated plots are saved as PDF files in the working directory.
#'
#' @param seurat_objects A named list of Seurat objects.
#' @param genes A character vector of gene names for which to generate feature plots.
#' @param filename The filename (without extension) to use for the output PDF files.
#'   Defaults to today's date.
#' @param separate Logical, whether to save a separate PDF for each Seurat object.
#'   If FALSE, all objects are combined into a single PDF. Default is TRUE.
#'
#' @return This function does not return a value. It generates PDF files with feature plots.
#'
#' @importFrom Seurat DefaultAssay FeaturePlot FetchData
#' @importFrom ggplot2 scale_color_gradientn ggtitle
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices pdf dev.off
#'
#' @export
generate_scaled_FPs <- function(seurat_objects, genes,
                                filename = as.character(Sys.Date()),
                                separate = TRUE) {

  if (is.null(names(seurat_objects))) {
    stop("seurat_objects must be a named list.")
  }

  # Find the maximum expression value for each gene across all objects
  max_values <- sapply(genes, function(gene) {
    vals <- sapply(seurat_objects, function(seurat_obj) {
      DefaultAssay(seurat_obj) <- "RNA"
      if (gene %in% rownames(seurat_obj)) {
        max(FetchData(seurat_obj, vars = gene), na.rm = TRUE)
      } else {
        NA_real_
      }
    })
    if (all(is.na(vals))) {
      warning("Gene '", gene, "' not found in any Seurat object and will be skipped.")
      NA_real_
    } else {
      max(vals, na.rm = TRUE)
    }
  })

  # Compute grid layout from number of genes
  n <- length(unique(genes))
  setncol <- ceiling(sqrt(n))
  setnrow <- ceiling(n / setncol)

  CombinePlots <- list()

  for (seurat_obj_name in names(seurat_objects)) {
    seurat_obj <- seurat_objects[[seurat_obj_name]]
    DefaultAssay(seurat_obj) <- "RNA"

    p <- FeaturePlot(seurat_obj, ncol = setncol, features = genes, combine = TRUE)

    for (i in seq_along(genes)) {
      gene <- genes[[i]]
      if (!is.na(max_values[[gene]])) {
        p[[i]] <- p[[i]] + scale_color_gradientn(
          colours = c("lightgrey", "navy"),
          limits = c(0, max_values[[gene]])
        )
      }
      if (i == 1) {
        p[[i]] <- p[[i]] + ggtitle(paste(seurat_obj_name, genes[1]))
      }
    }

    if (separate) {
      pdf_filename <- paste0(filename, "_", seurat_obj_name, "_FPs.pdf")
      pdf(file = pdf_filename, width = 3.75 * setncol, height = 3 * setnrow)
      print(p)
      dev.off()
      message("Saved: ", pdf_filename)
    } else {
      CombinePlots[[seurat_obj_name]] <- p
    }
  }

  if (!separate) {
    pdf_filename <- paste0(filename, "_FPs.pdf")
    pdf(file = pdf_filename, width = 3.75 * setncol, height = 3 * setnrow)
    print(CombinePlots)
    dev.off()
    message("Saved: ", pdf_filename)
  }
}
