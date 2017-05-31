#' Remove samples without annotation data from counts object
#'
#' This function filters a counts object, removing samples that are not found in a corresponding
#' sample annotation object.
#' @param counts a matrix or data frame of gene expression counts. Must have samples in columns and genes in rows, with colnames as sample identifiers.
#' @param design a data frame of sample information. At minimum, must contain a column corresponding to sample identifiers matching column names of \code{counts}.
#' @param libID_col numeric index or character name of column in \code{design} containing sample identifers matching column names of \code{counts}.
#' @details Samples are removed from \code{counts} if their column name does not match an element in column \code{libID_col} of \code{design}.
#' @export
#' @return a matrix or data frame with the same number of rows as \code{counts}
#' @usage \code{design_filter_counts(counts, design, libID_col)}
design_filter_counts <- function(counts, design, libID_col) {
  keepCols <- colnames(counts) %in% design[,libID_col]
  counts <- counts[,keepCols]
}
