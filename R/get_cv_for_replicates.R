#' calculate coefficient of variation (CV) values across pseudo-replicates in single-cell data
#'
#' @param counts : counts matrix from single-cell experiments, genes by cells
#' @param meta : metadata of cells, cells by features, must have sample and cell types coulumns
#' @param q : number of replicates

#'
#' @return : a list of CV across pseudo-replicates expression for each gene
#' @export

#' @examples
#' library(VICE)
#' data(cmat)
#' data(cmeta)
#' cvlist <- get_cv_for_replicates(cmeta, cmat, 3)
## Caculating CV for each gene across replicates
get_cv_for_replicates <- function(meta, counts, q) {
  stable <- table(meta$sample, meta$celltype)

  cvlist <- list()
  for (i in 1:nrow(stable)) {
    cvforcell <- list()
    for (n in 1:ncol(stable)) {
      sampleid <- rownames(stable)[[i]]
      cellid <- colnames(stable)[[n]]
      num <- as.numeric(stable[i, n])
      print(paste0("Sample: ", sampleid))
      print(paste0("Cell type: ", cellid))

      rep <- list()
      cvtmp <- list()
      if (num < 3) {
        rep <- cvtmp <- "NA"
      } else {
        rawcounts <- counts[, meta$sample == sampleid & meta$celltype == cellid]

        if (num/q > 1000) {
          cnum <- c(1, 3, 5, 10, 30, 50, seq(100, 900, 100), seq(1000, num / q, 1000))
          for (j in 1:length(cnum)) {
            print(cnum[[j]])
            rep[[j]] <- build_rep(rawcounts, cnum[[j]], q)
            exp0 <- apply(rep[[j]], 1, sum)
            rep[[j]] <- rep[[j]][exp0 != 0, ]
            cvtmp[[j]] <- apply(rep[[j]], 1, get_cv)
            names(cvtmp[[j]]) <- rownames(rep[[j]])
          }
        }

        if (num/q > 100 & (num/q < 1000 | num/q == 1000)) {
          cnum <- c(1, 3, 5, 10, 30, 50, seq(100, num / q, 100))
          for (j in 1:length(cnum)) {
            print(cnum[[j]])
            rep[[j]] <- build_rep(rawcounts, cnum[[j]], q)
            exp0 <- apply(rep[[j]], 1, sum)
            rep[[j]] <- rep[[j]][exp0 != 0, ]
            cvtmp[[j]] <- apply(rep[[j]], 1, get_cv)
            names(cvtmp[[j]]) <- rownames(rep[[j]])
          }
        }


        if (num/q > 10 & (num/q < 100 | num/q == 100)) {
          cnum <- c(1, 3, 5, seq(10, num / q, 10))
          for (j in 1:length(cnum)) {
            print(cnum[[j]])
            rep[[j]] <- build_rep(rawcounts, cnum[[j]], q)
            exp0 <- apply(rep[[j]], 1, sum)
            rep[[j]] <- rep[[j]][exp0 != 0, ]
            cvtmp[[j]] <- apply(rep[[j]], 1, get_cv)
            names(cvtmp[[j]]) <- rownames(rep[[j]])
          }
        }

        if (num/q < 10 | num/q == 10) {
          cnum <- seq(1, num / q, 1)
          for (j in 1:length(cnum)) {
            print(cnum[[j]])
            rep[[j]] <- build_rep(rawcounts, cnum[[j]], q)
            exp0 <- apply(rep[[j]], 1, sum)
            rep[[j]] <- rep[[j]][exp0 != 0, ]
            cvtmp[[j]] <- apply(rep[[j]], 1, get_cv)
            names(cvtmp[[j]]) <- rownames(rep[[j]])
          }
        }
        names(cvtmp) <- cnum
      }

      cvforcell[[n]] <- cvtmp
    }
    names(cvforcell) <- colnames(stable)
    cvlist[[i]] <- cvforcell
  }
  names(cvlist) <- rownames(stable)
  return(cvlist)
}
