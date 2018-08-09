#' Filter, normalize, and/or convert counts object for single-cell RNA-seq data
#'
#' This function performs a number of operations on a counts object. It optionally filters out low-count genes. It
#' optionally normalizes, log2-transforms, and/or transposes counts.
#' @param counts a matrix or data frame of gene expression counts. Should have samples in columns and genes in rows.
#' @param min_count numeric, the minimum count for a library to be considered passing threshold for a given gene. Setting to NULL avoids filtering by minimum counts.
#' @param min_cpm numeric, the minimum counts per million for a library to be considered passing threshold for a given gene. Setting to NULL avoids filtering by minimum counts.
#' @param min_libs_perc numeric, the minimum percentage of libraries that must meet \code{min_count} or \code{min_cpm} threshold for a gene to be retained.
#' @param normalize logical, whether to normalize counts.
#' @param norm_method character, the method by which to normalize counts. Used only if \code{normalize} is TRUE. Options are "TMM", which uses \code{edgeR::calcNormFactors}; "deconvolution", which uses cell deconvolution as implemented in the package \code{scran}; or "lib_size", simple library size normalization using cpm. Unique partial matches are accepted.
#' @param log2_transform logical, whether to log2 transform the counts. Defaults to \code{FALSE}.
#' @param transpose logical, whether to transpose the matrix or data frame of counts.
#' @param ... (optional) parameters passed to normalization functions.
#' @details This function (optionally) utilizes \code{min_filter_counts}
#'  to filter the counts object. It then (optionally) normalizes the counts, using
#'  one of several methods applicable to single-cell RNA-seq data. It then (optionally)
#'  log2-transforms the counts and/or transposes the counts object. Finally, it returns
#'  the counts either as a data frame, or as a \code{DGEList} object.
#' @export
#' @return a data frame containing the processed counts.
#' @usage \code{
#' calc_norm_counts_scRNAseq(
#'   counts,
#'   min_count=NULL, min_cpm=NULL, min_libs_perc=0.15,
#'   normalize=TRUE, norm_method="TMM",
#'   log2_transform=FALSE, transpose=FALSE, return_DGEcounts=FALSE,
#'   ...)}
calc_norm_counts_scRNAseq <-
  function(
    counts,
    min_count=NULL, min_cpm=NULL, min_libs_perc=0.15,
    normalize=TRUE, norm_method="TMM",
    log2_transform=FALSE, transpose=FALSE,
    ...) {
    
    norm_method <- tolower(norm_method)
    norm_method <-
      match.arg(norm_method, choices=c("tmm", "deconvolution", "lib_size"))
    # filter to keep genes that have minimum counts (or cpm) in minimum percent of libraries
    if (!is.null(min_cpm) | !is.null(min_count))
      counts <- min_filter_counts(counts, min_count=min_count, min_cpm=min_cpm, min_libs_perc)
    
    # normalize using TMM
    if (normalize & (norm_method %in% "tmm")) {
      DGEcounts <-
        edgeR::calcNormFactors(
          edgeR::DGEList(counts),
          method="TMM", ...)
      normCounts <- edgeR::cpm(DGEcounts, normalized.lib.sizes=TRUE)
    }
    
    # normalize using the cell deconvolution algorithm
    if (normalize & (norm_method == "deconvolution")) {
      if (!requireNamespace("scran", quietly = TRUE))
        stop("Package \"scran\" needed for deconvolution normalization. ",
             "Please install it or select a different normalization method.",
             call. = FALSE)
      decon_norm_factors <-
        scran::computeSumFactors(as.matrix(counts), ...)
      normCounts <- as.data.frame(t(t(counts)/decon_norm_factors))
    }
    
    # normalize using library size
    if (normalize & (norm_method == "lib_size")) {
      normCounts <- edgeR::cpm(counts, ...)
    }
    
    if (normalize) {
      return(normCounts)
    } else return(counts)
  }
