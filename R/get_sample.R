#' select sample with the largest number of cell
#'
#' @param meta : metadata of cells
#'
#' @return a dataframe of the samples with largest number of cells for each cell type
#' @export
#'
#' @examples
#' library(VICE)
#' data(cmat)
#' data(cmeta)
#' a <- get_sample(cmeta)
get_sample <- function(meta) {
  sample_cell <- as.data.frame(rbind(table(meta$sample, meta$celltype))) # sample-celltype-table
  sid <- apply(sample_cell, 2, function(x) {
    rownames(sample_cell)[which.max(x)]
  })
  snum <- apply(sample_cell, 2, function(x) {
    x[which.max(x)]
  })
  return(data.frame(cbind(sid, snum)))
}
