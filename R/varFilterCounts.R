#' Filter low-variance genes
#'
#' This function filters a counts object, keeping the \code{n_keep} genes with the highest variance.
#' @param counts a matrix or data frame of gene expression counts.
#' @param n_keep integer, the maximum number of genes to retain
#' @param genes character string, the dimension of \code{counts} containing the genes. Accepted values are "rows", "columns", and partial matches. Defaults to "rows"
#' @param log2_transform boolean, whether to log2 transform the counts prior to variance filtering. If the counts have not already been transformed, should be set to \code{TRUE}. Defaults to \code{FALSE}.
#' @details Variance is calculated for all genes. Genes are then ranked by this variance, and the \code{n_keep} genes with highest variance are retained.
#' @export
#' @return a matrix or data frame with the same number of samples as \code{counts}, but possibly with fewer genes.
#' @usage \code{
#' varFilterCounts(counts, counts, n_keep=8000, genes="rows",
#'                 log2_transform=FALSE
#'                 )}
varFilterCounts <- function(counts, n_keep=8000, genes="rows",
                            log2_transform=FALSE) {
  genes <- match.arg(genes, choices=c("rows", "columns"))
  if (log2_transform) counts <- log2(counts+1)
  if (genes=="columns") counts <- t(counts)
  if (nrow(counts) <= n_keep) { # if there are not more than n_keep genes in counts, return it without trimming
    if (genes=="columns") counts <- t(counts)
    return(counts)
  }
  var.counts <- apply(counts, MARGIN=1, FUN=var)
  keepGenes.tmp <- order(var.counts, decreasing=TRUE)[1:n_keep]
  counts.variable <- counts[keepGenes.tmp,]
  counts.variable <- counts.variable[order(rownames(counts.variable)),]
  if (genes=="columns") counts.variable <- t(counts.variable)
  return(counts.variable)
}
