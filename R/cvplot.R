## Plotting CV vs replicate size
# cvlist: list of CV values from get_cv_for_replicates function
# sampleid: specific sample to plot
# cellid: spcific cell type to plot
#' plot CV vs replicate size
#'
#' @param cvlist : a list of CV values calculated from get_cv_for_replicates
#' @param sampleid : ID of a sample for plot
#' @param cellid : ID of a cell type for plot
#'
#' @return a boxplot of CV in sample
#' @export
#'
#' @examples
#' library(VICE)
#' data(cmat)
#' data(cmeta)
#' cvlist <- get_cv_for_replicates(cmeta, cmat, 3)
#' cvplot(cvlist, "s1", "c1")
cvplot <- function(cvlist, sampleid, cellid) {
  id1 <- which(names(cvlist) == sampleid)
  id2 <- which(names(cvlist[[id1]]) == cellid)
  plotdata <- cvlist[[id1]][[id2]]
  if (length(plotdata) == 1 & any(is.na(plotdata))) {
    print(paste0("No ", cellid, " in ", sampleid))
  } else {
    par(mar = c(4, 4, 4, 4))
    boxplot(plotdata, las = 2, xlab = "Number of cells in replicate", ylab = "CV", main = paste0("Sample: ", sampleid, " , ", "Cell type: ", cellid))
    abline(h = 0.1, col = "red")
    abline(h = 0.2, col = "blue")
  }
}
