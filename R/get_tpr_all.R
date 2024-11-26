#' Estinate true positive rate of DE genes based on simulated data
#'
#' @param ngene : number of genes
#' @param ncell : number of cells per sample
#' @param nsample : number of samples in each DE group
#' @param de.prob : DE probability
#' @param de.facloc : Location (meanlog) parameter for the differential expression factor log-normal distribution. Changing these parameters can result in more or less extreme differences between groups.
#' @param de.facScale : Scale (sdlog) parameter for the differential expression factor log-normal distribution. Changing these parameters can result in more or less extreme differences between groups.

#' @param nperm : iteration times
#'
#' @return Overall TPR of DE analysis
#' @export
#'
#' @examples
#' library(VICE)
#' test <- get_tpr_all(ngene = 2000, ncell = 100, nsample = 3, de.prob = 0.3, de.facScale = 0.1, de.facloc = 0.1, nperm = 10)
get_tpr_all <- function(ngene = 2000, nsample = 3, ncell = 100, de.prob = 0.3, de.facloc = 0.1, de.facScale = 0.1, nperm = 10) {
  library(splatter)
  library(scater)
  library(edgeR)
  library(scuttle)
  library(ggplot2)
  library(reshape2)

  de <- function(data, group) {
    dgList <- DGEList(counts = data, group = factor(group))
    dgList <- calcNormFactors(dgList, method = "TMM")
    designMat <- model.matrix(~group)
    dgList <- estimateDisp(dgList, design = designMat)
    fit <- glmFit(dgList, designMat)
    lrt <- glmLRT(fit, coef = 2)
    return(topTags(lrt, n = nrow(dgList))$table)
  }

  getcv <- function(x) {
    sd(x) / mean(x)
  }
  norm <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  scell <- nsample * ncell * 2

  TPR <- list()
  TPR_l <- list()
  TPR_h <- list()
  for (z in 1:nperm) {
    sim.groups <- splatSimulate(
      method = "groups",
      nGenes = ngene,
      batchCells = scell,
      de.prob = de.prob,
      de.facLoc = de.facloc,
      de.facScale = de.facScale,
      group.prob = c(0.5, 0.5), verbose = F
    )

    df <- counts(sim.groups)
    gene <- rowData(sim.groups)
    cell <- colData(sim.groups)


    group1 <- subset(cell, Group == "Group1")$Cell
    group2 <- subset(cell, Group == "Group2")$Cell
    all_groups <- list(group1, group2)
    subgroups <- lapply(all_groups, function(group) {
      split(group, cut(seq_along(group), breaks = nsample, labels = FALSE))
    })

    subdf <- list()
    for (i in 1:length(subgroups)) {
      dflist <- lapply(subgroups[[i]], function(x) {
        df[, match(x, colnames(df))]
      })
      subdf[[i]] <- Reduce(cbind, lapply(dflist, function(x) {
        apply(x, 1, sum)
      }))
    }
    psdf <- Reduce(cbind, subdf)
    colnames(psdf) <- c(paste0("g1_", 1:nsample), paste0("g2_", 1:nsample))
    de1 <- de(psdf, c(rep("g1", nsample), rep("g2", nsample)))


    de1 <- de1[match(gene$Gene, rownames(de1)), ]

    cv <- apply(psdf[, 1:3], 1, getcv)
    de1 <- de1[!is.na(cv), ]
    gene <- gene[!is.na(cv), ]
    gene$cv <- cv[!is.na(cv)]

    snr <- norm(abs(de1$logFC)) / norm(gene$cv)
    stat <- cbind(gene, de1)
    stat$snr <- snr
    t1 <- subset(stat, DEFacGroup1 != 1 & DEFacGroup2 != 1)$Gene
    t2 <- subset(stat, FDR < 0.05)$Gene
    TPR[[z]] <- length(intersect(t1, t2)) / length(t1)

    print(z)
  }
  output <- unlist(TPR)
  return(output)
  print(paste("The averaged TPR is", mean(output)))
}
