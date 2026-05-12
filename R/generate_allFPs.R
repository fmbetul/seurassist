#' Generate Feature Plots for given genes
#'
#' This function generates feature plots for each gene provided in \code{gene_list}
#' in a Seurat object, arranges them, and saves the result as a PDF file.
#'
#' @param SeuFile A Seurat object that contains the gene expression data.
#' @param gene_list A named list of character vectors of gene names.
#' @param assay The assay to use. Default is "RNA".
#' @param reduction The dimensional reduction to use for plotting. Default is "umap".
#'
#' @return Saves a PDF file to the output directory. Returns invisibly.
#'
#' @importFrom Seurat DefaultAssay FeaturePlot
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 labs
#' @importFrom grDevices pdf dev.off
#'
#' @export
generate_AllFPs <- function(SeuFile, gene_list, assay = "RNA", reduction = "umap") {

  Filename <- as.character(substitute(SeuFile))
  DefaultAssay(SeuFile) <- assay

  # Filter each group to genes present in the object
  gene_list_2 <- lapply(gene_list, function(genes) {
    genes[genes %in% rownames(SeuFile)]
  })

  FPList <- list()

  for (d in names(gene_list_2)) {
    Genes <- gene_list_2[[d]]

    if (length(Genes) == 0) {
      message("Skipping '", d, "': no genes found in object.")
      next
    }

    FPSinglePage <- list()
    FPSinglePage[[1]] <- FeaturePlot(SeuFile, Genes[1], reduction = reduction) +
      labs(title = paste(d, Genes[1]))

    if (length(Genes) > 1) {
      for (p in 2:length(Genes)) {
        FPSinglePage[[p]] <- FeaturePlot(SeuFile, Genes[p], reduction = reduction)
      }
    }

    FPList[[d]] <- ggarrange(plotlist = FPSinglePage, ncol = 5, nrow = 4)
  }

  if (!file.exists("output")) dir.create("output")

  out_path <- paste0("output/", Filename, "_AllFPs.pdf")
  pdf(out_path, width = 25, height = 22)
  print(FPList)
  dev.off()
  message("Saved: ", out_path)

  invisible(NULL)
}
