#' Calculate counts and percentages of a metadata column in a Seurat object
#'
#' @description
#' This function takes a Seurat object and a specified metadata column,
#' then returns a data frame showing the count and percentage of each unique value
#' in that column. It is useful for summarising cell type distributions or other
#' categorical metadata within single-cell datasets.
#'
#' @param object A Seurat object.
#' @param metadata A character string specifying the column name in the metadata
#'   to summarise (e.g., "celltype", "sample", or "cluster").
#' @param colname A character string giving the desired name for the category column
#'   in the output data frame. Defaults to \code{"CellType"}.
#' @param sort Logical, whether to sort the resulting data frame by \code{Count}
#'   in decreasing order. Defaults to \code{TRUE}.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{<colname>}{Unique values from the specified metadata column.}
#'   \item{Count}{Number of cells in each category.}
#'   \item{Percent}{Percentage of cells in each category (rounded to one decimal place).}
#' }
#'
#' @examples
#' \dontrun{
#' percent_df <- count_percent_by(seurat_object, metadata = "celltype", colname = "CellType")
#' head(percent_df)
#'
#' # Disable sorting
#' percent_df_unsorted <- count_percent_by(seurat_object, metadata = "sample", sort = FALSE)
#' }
#'
#' @importFrom dplyr rename mutate arrange desc
#' @importFrom magrittr %>%
#' @importFrom rlang := .data
#'
#' @export
count_percent_by <- function(object, metadata, colname = "CellType", sort = TRUE) {

  if (!metadata %in% colnames(object[[]])) {
    stop("Column '", metadata, "' not found in object metadata.")
  }

  df <- as.data.frame(table(object[[]][[metadata]])) %>%
    dplyr::rename(!!colname := "Var1", Count = "Freq") %>%
    dplyr::mutate(Percent = round(100 * .data$Count / sum(.data$Count), 1))

  if (sort) {
    df <- df %>% dplyr::arrange(dplyr::desc(.data$Count))
  }

  return(df)
}
