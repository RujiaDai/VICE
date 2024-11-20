
#' Estinate true positive rate of DE genes based on simulated data
#'
#' @param ngene : number of genes
#' @param ncell : number of cells per sample
#' @param nsample : number of samples in each DE group
#' @param de.prob : DE probability
#' @param de.facloc : Differential expression factors are produced from a log-normal distribution. Changing these parameters can result in more or less extreme differences between groups.
#' @param nperm : iteration times
#'
#' @return TPR of DE analysis
#' @export
#'
#' @examples
#' library(VICE)
#' test<-get_tpr(ngene=1000,ncell=100,nsample=3,de.prob=0.2,de.facloc=0.05,nperm=5)

get_tpr<-function(ngene,nsample,ncell,de.prob,de.facloc,nperm){
library(splatter)
library(scater)
library(edgeR)
library(scuttle)
library(ggplot2)
library(reshape2)

  de<-function(data,group){
    dgList <- DGEList(counts=data,group=factor(group))
    dgList <- calcNormFactors(dgList, method="TMM")
    designMat <- model.matrix(~group)
    dgList <- estimateDisp(dgList, design=designMat)
    fit <- glmFit(dgList, designMat)
    lrt <- glmLRT(fit, coef=2)
    return(topTags(lrt,n=nrow(dgList))$table)}

  getcv<-function(x){sd(x)/mean(x)}
  norm<-function(x){(x-min(x))/(max(x)-min(x))}
  scell<-nsample*ncell*2

  TPR<-list()
  TPR_l<-list()
  TPR_h<-list()
  for(z in 1:nperm){
    sim.groups = splatSimulate(
        method='groups',
        nGenes = ngene,
        batchCells = scell,
        de.prob = de.prob,
        de.facLoc =de.facloc,#=0.1 for large effect size
        de.facScale=0.001,
        group.prob = c(0.5,0.5), verbose = F
      )

      df<-counts(sim.groups)
      gene<-rowData(sim.groups)
      cell<-colData(sim.groups)


      group1<-subset(cell,Group=='Group1')$Cell
      group2<-subset(cell,Group=='Group2')$Cell
      all_groups <- list(group1, group2)
      subgroups <- lapply(all_groups, function(group) {
        split(group, cut(seq_along(group), breaks = nsample, labels = FALSE))
      })

      subdf<-list()
      for(i in 1:length(subgroups)){
        dflist<-lapply(subgroups[[i]],function(x){df[,match(x,colnames(df))]})
        subdf[[i]]<-Reduce(cbind,lapply(dflist,function(x){apply(x,1,sum)}))
      }
      psdf<-Reduce(cbind,subdf)
      colnames(psdf)<-c(paste0("g1_",1:nsample),paste0("g2_",1:nsample))
      de1<-de(psdf,c(rep("g1",nsample),rep('g2',nsample)))


      de1<-de1[match(gene$Gene,rownames(de1)),]

      cv<-apply(psdf[,1:nsample],1,getcv)
      de1<-de1[!is.na(cv),]
      gene<-gene[!is.na(cv),]
      gene$cv<-cv[!is.na(cv)]

      snr<-norm(abs(de1$logFC))/norm(gene$cv)
      stat<-cbind(gene,de1)
      stat$snr<-snr


      t1<-subset(stat,DEFacGroup1!=1&DEFacGroup2!=1&(snr==2|snr>2)&snr<4)$Gene
      t2<-subset(stat,FDR<0.05&(snr==2|snr>2)&snr<4)$Gene
      TPR[[z]]<-length(intersect(t1,t2))/length(t1)

      t3<-subset(stat,DEFacGroup1!=1&DEFacGroup2!=1&(snr>4|snr==4))$Gene
      t4<-subset(stat,FDR<0.05&(snr>4|snr==4))$Gene
      TPR_h[[z]]<-length(intersect(t4,t3))/length(t3)

      t5<-subset(stat,DEFacGroup1!=1&DEFacGroup2!=1&snr>1&snr<2)$Gene
      t6<-subset(stat,FDR<0.05&snr>1&snr<2)$Gene
      TPR_l[[z]]<-length(intersect(t5,t6))/length(t5)

    print(z)
  }
  output<-data.frame(high=unlist(TPR_h),medium=unlist(TPR),low=unlist(TPR_l))
  meanfile<-data.frame(cbind(apply(output,2,function(x){mean(x,na.rm=T)}),apply(output,2,function(x){sd(x,na.rm=T)})))
  colnames(meanfile)<-c("mean","sd")
  meanfile$variable<-rownames(meanfile)
  meanfile$variable<-factor(meanfile$variable,levels=c("high","medium","low"))

  g1<-ggplot(meanfile, aes(x=variable, y=mean, group=variable, color=variable)) +
    geom_point(position=position_dodge(0.05),size=4)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05))+
    ggpubr::theme_pubclean()+ylim(0,1)+
    labs(x='SNR',y='TPR')+
    scale_color_manual(values = c("#de2d26", "#fc9272","#fee0d2"))+
	theme(legend.position='none')
  print(g1)
  return(output)
}



