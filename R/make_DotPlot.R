# make_DotPlot(): DotPlot wrapper that groups genes into named panels with
# dividers, dedupes genes appearing in several panels, and reverses the group.by
# levels so the first level plots at the top. Adapted from the 28-novo-in-vivo
# pipeline.

#' Build a grouped DotPlot from a Seurat object
#'
#' Wraps [Seurat::DotPlot()], splitting genes into named panels, dropping genes
#' absent from the assay, and optionally de-duplicating genes across panels.
#' The caller's default assay and identities are restored on exit.
#'
#' @param seurat_obj    Seurat object.
#' @param gene_list     Named list of character vectors, one per panel.
#' @param group_by      Metadata column to set as Idents(); NA uses the
#'                      object's current Idents().
#' @param assay         Assay to read expression from.
#' @param ensure_unique  Drop a gene from every panel after its first
#'                      occurrence, so DotPlot does not draw it twice.
#' @param reverse_levels Reverse group_by's factor levels so the first level
#'                      plots at the top of the y axis.
#' @param palette       Colours for [ggplot2::scale_colour_gradientn()].
#' @param base_size     Base font size.
#' @param base_family   Font family (default `"sans"`).
#' @param angle_x       Rotation of the x-axis (gene) labels.
#' @return list(plot = ggplot, groups_kept = named int vector of genes kept
#'   per panel after filtering to those present in the assay).
#' @importFrom Seurat DefaultAssay Idents DotPlot
#' @examples
#' \dontrun{
#' make_DotPlot(obj, list(Neurons = c("Snap25", "Rbfox3"), Glia = c("Gfap")),
#'              group_by = "seurat_clusters")$plot
#' }
#' @export
make_DotPlot <- function(
    seurat_obj,
    gene_list,
    group_by       = NA_character_,
    assay          = "RNA",
    ensure_unique  = TRUE,
    reverse_levels = TRUE,
    palette        = c("#2E7094", "#FFFAE8", "#B00012"),
    base_size      = 12,
    base_family    = "sans",
    angle_x        = 45
) {

  stopifnot(is.list(gene_list), length(gene_list) > 0)

  if (!assay %in% names(seurat_obj@assays))
    stop(sprintf("Assay '%s' not found in this Seurat object. Available: %s",
                 assay, paste(names(seurat_obj@assays), collapse = ", ")))

  # Preserve caller's assay/idents.
  old_assay  <- DefaultAssay(seurat_obj)
  old_idents <- Idents(seurat_obj)
  on.exit({
    DefaultAssay(seurat_obj) <- old_assay
    Idents(seurat_obj)       <- old_idents
  }, add = TRUE)

  DefaultAssay(seurat_obj) <- assay

  if (!is.na(group_by)) {
    if (!group_by %in% colnames(seurat_obj@meta.data))
      stop(sprintf("Metadata column '%s' not found. Available columns: %s",
                   group_by, paste(colnames(seurat_obj@meta.data), collapse = ", ")))
    vals <- seurat_obj[[group_by, drop = TRUE]]
    lvls <- if (is.factor(vals)) levels(vals) else
      tryCatch(as.character(sort(as.numeric(unique(as.character(vals))))),
               warning = function(w) sort(unique(as.character(vals))))
    Idents(seurat_obj) <- factor(as.character(vals), levels = lvls)
  }

  if (reverse_levels && is.factor(Idents(seurat_obj)))
    levels(seurat_obj) <- rev(levels(seurat_obj))

  genes_present  <- rownames(seurat_obj[[assay]])
  used           <- character(0)
  gene_list_filt <- lapply(gene_list, function(g) {
    g_in <- intersect(unique(g), genes_present)
    if (ensure_unique) {
      g_in <- setdiff(g_in, used)
      used <<- c(used, g_in)
    }
    g_in
  })
  gene_list_filt <- gene_list_filt[vapply(gene_list_filt, length, 1L) > 0]

  if (length(gene_list_filt) == 0)
    stop("No genes remain after filtering. Check that gene names match the assay's feature names.")

  message(sprintf("[make_DotPlot] Plotting %d gene groups | %d total genes",
                  length(gene_list_filt),
                  sum(vapply(gene_list_filt, length, 1L))))

  p <- Seurat::DotPlot(seurat_obj, features = gene_list_filt, scale = TRUE) +
    ggplot2::scale_colour_gradientn(colours = palette) +
    ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = angle_x, hjust = 1, vjust = 1,
                                           color = "gray10", face = "italic",
                                           size  = base_size),
      axis.text.y  = ggplot2::element_text(color = "black", size = base_size,
                                           face  = "bold"),
      legend.title = ggplot2::element_text(size = base_size - 2),
      legend.text  = ggplot2::element_text(size = base_size - 2),
      title        = ggplot2::element_text(size = base_size),
      axis.title   = ggplot2::element_text(size = base_size - 1, face = "plain")
    ) +
    ggplot2::xlab("Genes") + ggplot2::ylab("Cluster")

  list(plot = p, groups_kept = vapply(gene_list_filt, length, 1L))
}
