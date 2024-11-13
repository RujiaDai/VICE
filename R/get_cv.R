#' Caculating coefficient of variation
#'
#' @param x a vector of gene expression
#'
#' @return CV value
#' @export
#'
#' @examples
#' library(VICE)
#' data(cmat)
#' data(cmeta)
#' get_cv(cmat[1, ]) # calculate CV for the first gene
get_cv <- function(x) {
  sd(x) / mean(x)
}
