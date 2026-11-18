#' Save QC file
#'
#' This function saves UMAP plot and quality control (QC) files, with options
#' to include QC plots and heatmaps, to the output folder.
#'
#' @param SeuFile Seurat object that contains the data.
#' @param version A string that provides the version number for the QC file. Default is "v1.0.0".
#' @param qc.plot Logical, if TRUE the function will generate QC plots for each cluster. Default is TRUE.
#' @param Heatmap Logical, if TRUE the function will generate a heatmap of top marker genes. Default is FALSE.
#' @param csv Logical, if TRUE the function will write CSV files for the top 10 markers and all markers
#'   per cluster. Requires Heatmap = TRUE (marker genes must be computed first). Default is FALSE.
#' @param pdf.height A numeric value specifying the height of the PDF. Default is 10.
#' @param pdf.width A numeric value specifying the width of the PDF. Default is 12.
#'
#' @return The function does not return anything. It saves the PDF and optionally CSV files in
#'   the "output" directory and prints plots to the console.
#'
#' @seealso \code{\link[Seurat]{DimPlot}}, \code{\link[Seurat]{VlnPlot}}, \code{\link[Seurat]{DoHeatmap}}
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.csv
#' @importFrom Seurat DimPlot PercentageFeatureSet VlnPlot FindAllMarkers DoHeatmap
#' @importFrom dplyr group_by top_n
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
save_qc_file <- function(SeuFile, version = "v1.0.0", qc.plot = TRUE, Heatmap = FALSE,
                         csv = FALSE, pdf.height = 10, pdf.width = 12) {

  Filename <- as.character(substitute(SeuFile))

  if (!file.exists("output")) dir.create("output")

  out_path <- paste0("output/", Filename, "_qc_file_", version, ".pdf")
  pdf(out_path, width = pdf.width, height = pdf.height)

  # UMAP
  DefaultAssay(SeuFile) <- "RNA"
  umap <- DimPlot(SeuFile, reduction = "umap", label = TRUE, pt.size = 0.5)
  print(umap)

  if (qc.plot) {
    SeuFile[["percent.mt"]] <- PercentageFeatureSet(SeuFile, pattern = "^MT-")
    qc <- VlnPlot(SeuFile, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(qc)
  }

  if (Heatmap) {
    message("Running FindAllMarkers - this may take a while...")
    markers <- FindAllMarkers(SeuFile, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    top10 <- markers %>%
      dplyr::group_by(.data$cluster) %>%
      dplyr::top_n(n = 10, wt = .data$avg_log2FC)

    hm <- DoHeatmap(SeuFile, features = top10$gene)
    print(hm)

    if (csv) {
      write.csv(top10,
                file = paste0("output/", Filename, "_MarkersPerCluster_top10_", version, ".csv"))
      write.csv(markers,
                file = paste0("output/", Filename, "_MarkersPerCluster_", version, ".csv"))
      message("Saved CSV files to output/")
    }
  } else if (csv) {
    warning("csv = TRUE has no effect when Heatmap = FALSE: markers are only computed when Heatmap = TRUE.")
  }

  dev.off()
  message("Saved: ", out_path)

  # Print to console for interactive use
  print(umap)
  if (qc.plot) print(qc)
  if (Heatmap) print(hm)
}
