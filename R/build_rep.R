#' construct pseudo-replicates in single-cell data
#'
#' @param counts : counts matrix from single-cell experiments, genes by cells
#' @param m : replicate size, number of cells in each replicate
#' @param q : number of replicates

#' @return a list of pseudo-replicates
#' @export
#'
#' @examples
#' library(VICE)
#' data(cmat)
#' data(cmeta)
#' psrep <- build_rep(cmat, 100, 3)
build_rep <- function(counts, m, q) {
  replicate_list <- c()
  index <- colnames(counts)
  k <- 1
  while (k <= q) {
    fz <- sample(index, m, replace = FALSE, prob = NULL)
    counts_fz <- counts[, fz]
    if (m >= 2) {
      replicate_i <- apply(counts_fz, 1, sum)
      replicate_list <- cbind(replicate_list, replicate_i)
      rownames(replicate_list) <- rownames(counts)
    } else {
      replicate_list <- cbind(replicate_list, counts_fz)
      rownames(replicate_list) <- rownames(counts)
    }
    for (i in 1:m) {
      index <- index[-which(index == fz[i])]
    }
    k <- k + 1
  }

  return(replicate_list)
}
