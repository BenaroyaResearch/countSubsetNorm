#' Trim excess factor levels from columns of a data frame
#'
#' This function goes through all columns of a data frame, and removes excess levels from columns of
#' class \code{factor}, while leaving other columns unchanged. It is intended to be used after
#' filtering a data frame, so that excess levels do not cause errors or confusion in downstream
#' processes.
#' @param data  a data frame, with 0 or more columns of class \code{factor}
#' @export
#' @return a data frame of the same dimensions as \code{data}
#' @usage \code{fix_factors(data)}
fix_factors <- function(data) {
  data <- as.data.frame(lapply(data, FUN=function(x) {if (is.factor(x)) factor(x) else x}))
}