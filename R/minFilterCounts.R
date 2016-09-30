#' Filter low-count genes
#'
#' This function filters a counts object, removing genes that do not meet a minimum count value in
#' a minimum percentage of libraries
#' @param counts a matrix or data frame of gene expression counts; should have samples in columns and genes in rows.
#' @param min_count numeric value, the minimum count value
#' @param min_libs_perc numeric value, the minimum percentage of libraries
#' @details Genes are removed from \code{counts} if they do not have at least \code{min_count} counts in at least \code{min_libs_perc} % of all libraries.
#'  Note that adding or removing libraries to/from the counts object may change the genes that meet these thresholds.
#' @export
#' @return a matrix or data frame with the same number of columns as \code{counts}
#' @usage \code{minFilterCounts(counts, min_count=1, min_libs_perc=0.15)}
minFilterCounts <- function(counts, min_count=1, min_libs_perc=0.15) {
  keepRows <- rowSums((counts) >= min_count) >= (min_libs_perc * ncol(counts))
  counts <- counts[keepRows,]
}