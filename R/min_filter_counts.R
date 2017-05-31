#' Filter low-count genes
#'
#' This function filters a counts object, removing genes that do not meet a minimum count value in
#' a minimum percentage of libraries
#' @param counts a matrix or data frame of gene expression counts; should have samples in columns and genes in rows.
#' @param min_count numeric value, the minimum count value. Either this or min_cpm should be specified, but not both.
#' @param min_cpm numeric values, the minimum counts per million value. Either this or min_cpm should be specified, but not both.
#' @param min_libs_perc numeric value, the minimum percentage of libraries with counts or cpm equal to or greater than the threshold.
#' @details Genes are removed from \code{counts} if they do not have at least \code{min_count} counts in at least \code{min_libs_perc} % of all libraries.
#'  Note that adding or removing libraries to/from the counts object may change the genes that meet these thresholds.
#' @export
#' @return a filtered matrix or data frame with the same number of columns as \code{counts}, and potentially fewer rows. Counts are retained in the units they are input (using min_cpm does not convert counts to cpm).
#' @usage \code{min_filter_counts(counts, min_count=NULL, min_cpm=NULL, min_libs_perc=0.15)}
min_filter_counts <-
  function(counts, min_count=NULL, min_cpm=NULL, min_libs_perc=0.15) {
  if (!is.null(min_count)) {
    if (!is.null(min_cpm)) {
      stop("Both minimum cpm and minimum count were specified; please choose one or the other.")
    } else {
      keepRows <- rowSums((counts) >= min_count) >= (min_libs_perc * ncol(counts))
    }
  } else if (!is.null(min_cpm)) {
    keepRows <- rowSums((edgeR::cpm(counts)) >= min_cpm) >= (min_libs_perc * ncol(counts))
  } else stop("Neither minimum cpm or minimum count were specified; please set min_count or min_cpm.")
  counts <- counts[keepRows,]
}
