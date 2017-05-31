#' Extract counts matrix from different types of expression objects
#'
#' This is a helper function to extract the counts matrix from different gene expression object classes. It is used
#' by many other functions to enable counts being input as matrices, data.frames, EList or DGEList objects. It
#' determines the class of the object and returns the counts as a data.frame (default) or matrix.
#' @param counts an object from which counts can be extracted. May be of class matrix, data.frame, EList, DGElist (others to be added).
#' @param return_class optional, the class for the object to be returned. Defaults to NULL, which causes the counts portion of the object to be returned as-is.
#' @export
#' @return A matrix or data frame.
#' @usage \code{
#' extract_counts(counts, return_class=NULL)}
extract_counts <-
  function(counts, return_class=NULL) {
    input_class <- class(counts)
    
    if ("EList" %in% input_class) {
      counts <- counts[["E"]]
    } else if ("DGEList" %in% input_class) {
      counts <- counts[["counts"]]
    } else if (("ExpressionSet" %in% input_class) | ("eSet" %in% input_class)) {
      if (is.environment(counts@assayData)) {
        counts <- counts@assayData[["exprs"]]
      } else if (is.matrix(counts@assayData)) {
        counts <- counts@assayData
      }
    }
    
    if (!is.null(return_class)) counts <- as(counts, return_class)
    
    counts
  }
