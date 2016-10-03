#' Filter, normalize, and/or convert counts object
#'
#' This function performs a number of operations on a counts object. It filters the counts to include
#' only samples matching a sample annotation object. It optionally filters out low-count genes. It
#' optionally normalizes, log2-transforms, and/or transposes counts. Finally, it can return either
#' a DGEList object or a data frame.
#' @param counts a matrix or data frame of gene expression counts. Should have sample in columns and genes in rows.
#' @param design a data frame of sample information. At minimum, must contain a column corresponding to sample identifiers matching column names of \code{counts}. Passed to \code{designFilterCounts}.
#' @param libID_col numeric index or character name of column in \code{design} containing sample identifers matching column names of \code{counts}. Passed to \code{designFilterCounts}.
#' @param min_count numeric, the minimum count for a library to be considered passing threshold for a given gene.
#' @param min_cpm numeric, the minimum counts per million for a library to be considered passing threshold for a given gene.
#' @param min_libs_perc numeric, the minimum percentage of libraries that must meet \code{min_count} or \code{min_cpm} threshold for a gene to be retained.
#' @param normalize logical, whether to normalize counts using \code{edgeR::calcNormFactors}.
#' @param norm_method character, the method by which to normalize counts; passed to \code{edgeR::calcNormFactors}. Used only if \code{normalize} is TRUE. Defaults to "TMM".
#' @param log2_transform logical, whether to log2 transform the counts. Defaults to \code{FALSE}.
#' @param transpose logical, whether to transpose the matrix or data frame of counts.
#' @param return_DGEcounts logical, whether to return return counts as a \code{DGEList} object. If FALSE, counts are returned as a data frame. Defaults to FALSE.
#' @param ... (optional) parameters passed to \code{DGEList}.
#' @details This function utilizes \code{designFilterCounts} and (optionally) \code{minFilterCounts}
#'  to filter the counts object. It then (optionally) normalizes the counts, using
#'  \code{edgeR::calcNormFactors}. It then (optionally) log2-transforms the counts and/or transposes
#'  the counts object. Finally, it returns the counts either as a data frame, or as a \code{DGEList}
#'  object.
#'  For WGCNA analyses, should use \code{log2_transform=TRUE, transpose=TRUE, return_DGECounts=FALSE}.
#'  For limma analyses, should use \code{log2_transform=FALSE, transpose=FALSE, return_DGECounts=TRUE}.
#' @export
#' @return a data frame or \code{DGEList} object containing the processed counts.
#' @usage \code{
#' calcNormCounts(counts, design, libID_col="lib.id",
#'                min_count=NULL, min_cpm=NULL, min_libs_perc=0.15,
#'                normalize=TRUE, norm_method="TMM",
#'                log2_transform=FALSE, transpose=FALSE, return_DGEcounts=FALSE,
#'                ...
#'                )}
calcNormCounts <- function(counts, design, libID_col="lib.id",
                           min_count=NULL, min_cpm=NULL, min_libs_perc=0.15,
                           normalize=TRUE, norm_method="TMM",
                           log2_transform=FALSE, transpose=FALSE, return_DGEcounts=FALSE,
                           ...) {
  # trim counts object to include only desired libraries
  counts <- designFilterCounts(counts, design, libID_col)
  
  # filter to keep genes that have minimum counts (or cpm) in minimum percent of libraries
  if (!is.null(min_cpm) | !is.null(min_count))
    counts <- minFilterCounts(counts, min_count=min_count, min_cpm=min_cpm, min_libs_perc)
  
  # generate DGEList object, and normalize counts
  if (return_DGEcounts | normalize) {
    DGEcounts <- edgeR::DGEList(counts,...)
    if (normalize) DGEcounts <- edgeR::calcNormFactors(DGEcounts, method=norm_method)
  }
  if (return_DGEcounts) return(DGEcounts)
  
  if (normalize) {
    normCounts <- as.data.frame(sweep(counts, 2, DGEcounts$samples$norm.factors, "/"))
  } else normCounts <- counts
  
  # log-transform and/or transpose counts
  if (log2_transform) normCounts <- log2(normCounts + 1)
  if (transpose) normCounts <- as.data.frame(t(normCounts))
  
  return(normCounts)
}
