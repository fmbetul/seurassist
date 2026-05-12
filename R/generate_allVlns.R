#' Generate Violin Plots for given genes across different clusters
#'
#' This function generates violin plots for each gene provided in \code{gene_list}
#' across different clusters in a Seurat object, arranges them, and saves the result
#' as a PDF file.
#'
#' @param SeuFile A Seurat object that contains the gene expression data.
#' @param gene_list A named list of character vectors of gene names.
#' @param assay The assay to use. Default is "RNA".
#'
#' @return Saves a PDF file to the output directory. Returns invisibly.
#'
#' @importFrom Seurat DefaultAssay VlnPlot Idents
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 coord_flip theme element_blank element_text labs
#' @importFrom grDevices pdf dev.off
#'
#' @export
generate_allVlns <- function(SeuFile, gene_list, assay = "RNA") {

  Filename <- as.character(substitute(SeuFile))
  DefaultAssay(SeuFile) <- assay

  # Filter each group to genes present in the object
  gene_list_2 <- lapply(gene_list, function(genes) {
    genes[genes %in% rownames(SeuFile)]
  })

  AllVlns <- list()

  for (d in names(gene_list_2)) {
    Genes <- gene_list_2[[d]]

    if (length(Genes) == 0) {
      message("Skipping '", d, "': no genes found in object.")
      next
    }

    StorePlots <- list()

    # First gene: carries the group label as title
    plotA <- VlnPlot(SeuFile, features = Genes[1], pt.size = 0, same.y.lims = FALSE)
    plotA <- plotA + coord_flip() +
      theme(
        axis.ticks.x  = element_blank(),
        axis.title.y   = element_blank(),
        axis.text.x    = element_blank(),
        axis.title.x   = element_blank(),
        axis.ticks.y   = element_blank(),
        legend.position = "none",
        plot.title     = element_text(size = 12)
      ) +
      labs(title = d, subtitle = Genes[1])
    StorePlots[[Genes[1]]] <- plotA

    # Remaining genes: hide y-axis labels to save space
    if (length(Genes) > 1) {
      for (x in Genes[2:length(Genes)]) {
        plotB <- VlnPlot(SeuFile, features = x, pt.size = 0, same.y.lims = FALSE)
        plotB <- plotB + coord_flip() +
          theme(
            axis.ticks.x   = element_blank(),
            axis.title.y    = element_blank(),
            axis.text.x     = element_blank(),
            axis.title.x    = element_blank(),
            axis.ticks.y    = element_blank(),
            axis.text.y     = element_blank(),
            legend.position = "none",
            plot.title      = element_text(size = 12)
          )
        StorePlots[[x]] <- plotB
      }
    }

    AllVlns[[d]] <- ggarrange(
      plotlist = StorePlots,
      widths   = c(1.4, rep(1, length(Genes) - 1)),
      ncol     = 20,
      nrow     = 1
    )
  }

  if (!file.exists("output")) dir.create("output")

  out_path <- paste0("output/", Filename, "_AllVlns.pdf")
  pdf(out_path, width = 40, height = length(levels(Idents(SeuFile))))
  print(AllVlns)
  dev.off()
  message("Saved: ", out_path)

  invisible(NULL)
}
