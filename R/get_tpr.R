
#' Estinate true positive rate of DE genes based on simulated data
#'
#' @param ngene : number of genes
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
#' test<-get_tpr(ngene=500,nsample=3,de.prob=0.2,de.facloc=0.05,nperm=5)

get_tpr<-function(ngene,nsample,de.prob,de.facloc,nperm){
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

  scell<-c(30,150,300,1500,3000,4500,6000,9000,12000,15000)*2

  TPR<-list()
  TPR_l<-list()
  TPR_h<-list()
  for(z in 1:nperm){

    tpr<-c()
    tpr_h<-c()
    tpr_l<-c()
    for(g in 1:length(scell)){
      sim.groups = splatSimulate(
        method='groups',
        nGenes = ngene,
        batchCells = scell[[g]],
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


      t1<-subset(stat,DEFacGroup1!=1&DEFacGroup2!=1)$Gene
      t2<-subset(stat,FDR<0.05)$Gene
      tpr[[g]]<-length(intersect(t1,t2))/length(t1)


      t1<-subset(stat,DEFacGroup1!=1&DEFacGroup2!=1&(snr==2|snr>2)&snr<4)$Gene
      t2<-subset(stat,FDR<0.05&(snr==2|snr>2)&snr<4)$Gene
      tpr[[g]]<-length(intersect(t1,t2))/length(t1)

      t3<-subset(stat,DEFacGroup1!=1&DEFacGroup2!=1&(snr>4|snr==4))$Gene
      t4<-subset(stat,FDR<0.05&(snr>4|snr==4))$Gene
      tpr_h[[g]]<-length(intersect(t4,t3))/length(t3)

      t5<-subset(stat,DEFacGroup1!=1&DEFacGroup2!=1&snr>1&snr<2)$Gene
      t6<-subset(stat,FDR<0.05&snr>1&snr<2)$Gene
      tpr_l[[g]]<-length(intersect(t5,t6))/length(t5)

    }

    TPR[[z]]<-tpr
    TPR_l[[z]]<-tpr_l
    TPR_h[[z]]<-tpr_h
    print(z)
  }

  output<-matrix(as.numeric(Reduce(rbind,TPR)),nrow=nrow(Reduce(rbind,TPR)),ncol=ncol(Reduce(rbind,TPR)))
  colnames(output)<-scell/6
  rownames(output)<-1:nperm


  output2<-matrix(as.numeric(Reduce(rbind,TPR_l)),nrow=nrow(Reduce(rbind,TPR_l)),ncol=ncol(Reduce(rbind,TPR_l)))
  colnames(output2)<-scell/6
  rownames(output2)<-1:nperm

  output3<-matrix(as.numeric(Reduce(rbind,TPR_h)),nrow=nrow(Reduce(rbind,TPR_h)),ncol=ncol(Reduce(rbind,TPR_h)))
  colnames(output3)<-scell/6
  rownames(output3)<-1:nperm


  meanfile<-data.frame(high=unlist(apply(output3,2,function(x){mean(x[!is.na(x)])})),medium=unlist(apply(output,2,function(x){mean(x[!is.na(x)])})),low=unlist(apply(output2,2,function(x){mean(x[!is.na(x)])})))
  sdfile<-data.frame(high=unlist(apply(output3,2,function(x){sd(x[!is.na(x)])})),medium=unlist(apply(output,2,function(x){sd(x[!is.na(x)])})),low=unlist(apply(output2,2,function(x){sd(x[!is.na(x)])})))
  meanfile$num<-rownames(meanfile)
  meanfile<-melt(meanfile,id='num')
  colnames(meanfile)[3]<-"mean"
  meanfile$sd<-melt(sdfile)$value
  meanfile$num<-factor(meanfile$num,levels=as.factor(scell/6))


  ggplot(meanfile, aes(x=num, y=mean, group=variable, color=variable)) +
    geom_point(position=position_dodge(0.05),size=4)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05))+
    ggpubr::theme_pubclean()+ylim(0,1)+
    labs(x='Cell number',y='TPR')+
    scale_color_manual(values = c("#de2d26", "#fc9272","#fee0d2"))
  meanfile
}



