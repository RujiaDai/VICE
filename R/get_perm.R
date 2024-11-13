#' permutation test for CV
#'
#' @param meta : metadata of cells
#' @param counts : counts matrix from single-cell experiments, genes by cells
#' @param cvlist :cvlist calculated from get_cv_for_replicates
#' @param q : number of replicates
#' @param k : permutation times
#'
#' @return a list of CVs
#' @export
#'
#' @examples
#' library(VICE)
#' data(cmat)
#' data(cmeta)
#' cvlist <- get_cv_for_replicates(cmeta, cmat, 3)
#' permcv <- get_perm(cmeta, cmat, cvlist, 3, 10)
get_perm <- function(meta, counts, cvlist, q, k) {
  stable <- table(meta$sample, meta$celltype)
  pval <- list()
  testcv <- list()
  for (j in 1:ncol(stable)) {
    sampleid <- rownames(stable)[which.max(stable[, j])]
    cellid <- colnames(stable)[[j]]
    num <- names(cvlist[[sampleid]][[cellid]])[[length(names(cvlist[[sampleid]][[cellid]]))]]
    num <- as.numeric(num)

    rawcounts <- counts[, meta$sample == sampleid & meta$celltype == cellid]
    cvtmp <- list()
    for (i in 1:k) {
      df <- build_rep(rawcounts, num, q)

      cvtmp[[i]] <- apply(df, 1, get_cv)
    }
    cvtmp2 <- Reduce(cbind, cvtmp)

    mincv <- cvlist[[sampleid]][[cellid]][[as.character(num)]]
    testcv[[j]] <- cvtmp2[match(names(mincv), rownames(cvtmp2)), ]
  }
  names(testcv) <- colnames(stable)
  testcv
  # p<-c()
  # for(n in 1:length(mincv)){
  # testcv<-cvtmp2[match(names(mincv)[[n]],rownames(cvtmp2)),]
  # upper<-mean(testcv)+3*sd(testcv)
  # lower<-mean(testcv)-3*sd(testcv)
  # p[[n]]<-ifelse(mincv[[n]]<upper&mincv[[n]]>lower,1,0)
  # }
  # pval[[j]]<-unlist(p)
  # names(pval[[j]])<-names(mincv)
  # }
  # names(pval)<-colnames(stable)
  # pval
}
