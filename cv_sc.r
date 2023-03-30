#VICE: variability in single-cell/neclei RNAseq data


##Building replicates
#n:number of replicates
#m:replicate size, number of cells in each replicate
build_rep<- function(counts,m,n){
replicate_list=c()
index=colnames(counts)
k=1  
while (k<=n){
	fz <- sample(index,m,replace = FALSE,prob = NULL)
	counts_fz=counts[,fz]
	#sample_meta_fz=sample_meta[fz,]
	if(m>=2){
		replicate_i=apply(counts_fz,1,sum)
		replicate_list=cbind(replicate_list,replicate_i)}
	else{replicate_list=cbind(replicate_list,counts_fz)}
	for (i in 1:m) {
		index <- index[-which(index==fz[i])]
	}
k=k+1}

return(replicate_list)
}

##Caculating coefficient of variation 
get_cv=function(x){
sd(x)/mean(x)
}

##Caculating CV for each gene across replicates
#meta:meta matrix, cells by features, must have sample and cell types coulumns
#counts:read count matrix, genes by cells
get_cv_for_replicates<-function(meta,counts){
stable<-table(meta$sample,meta$celltype)

cvlist<-list()
for(i in 1:nrow(stable)){
cvforcell<-list()
for(n in 1:ncol(stable)){

sampleid=rownames(stable)[[i]]
cellid<-colnames(stable)[[n]]
num=as.numeric(stable[i,n])
print(paste0('Sample: ',sampleid))
print(paste0('Cell type: ',cellid))

rep<-list()
cvtmp<-list()
if(num<3){rep<-cvtmp<-'NA'}
else{
rawcounts<-counts[,meta$sample==sampleid&meta$celltype==cellid]  

if(num>3000){
cnum<-c(1,3,5,10,30,50,seq(100,900,100),seq(1000,num/3,1000))
for(j in 1:length(cnum)){
print(cnum[[j]])
rep[[j]]<-build_rep(rawcounts,cnum[[j]],3)
exp0<-apply(rep[[j]],1,sum)
rep[[j]]<-rep[[j]][exp0!=0,]
cvtmp[[j]]<-apply(rep[[j]],1,get_cv)
}}

if(num>300&(num<3000|num==3000)){
cnum<-c(1,3,5,10,30,50,seq(100,num/3,100))
for(j in 1:length(cnum)){
print(cnum[[j]])
rep[[j]]<-build_rep(rawcounts,cnum[[j]],3)
exp0<-apply(rep[[j]],1,sum)
rep[[j]]<-rep[[j]][exp0!=0,]
cvtmp[[j]]<-apply(rep[[j]],1,get_cv)
}}


if(num>30&(num<300|num==300)){
cnum<-c(1,3,5,seq(10,num/3,10))
for(j in 1:length(cnum)){
print(cnum[[j]])
rep[[j]]<-build_rep(rawcounts,cnum[[j]],3)
exp0<-apply(rep[[j]],1,sum)
rep[[j]]<-rep[[j]][exp0!=0,]
cvtmp[[j]]<-apply(rep[[j]],1,get_cv)
}}

if(num<30|num==30){
cnum<-seq(1,num/3,1)
for(j in 1:length(cnum)){
print(cnum[[j]])
rep[[j]]<-build_rep(rawcounts,cnum[[j]],3)
exp0<-apply(rep[[j]],1,sum)
rep[[j]]<-rep[[j]][exp0!=0,]
cvtmp[[j]]<-apply(rep[[j]],1,get_cv)
}}
names(cvtmp)<-cnum
}

cvforcell[[n]]<-cvtmp
}
names(cvforcell)<-colnames(stable)
cvlist[[i]]<-cvforcell
}
names(cvlist)<-rownames(stable)
return(cvlist)
}

##Plotting CV vs replicate size
#cvlist: list of CV values from get_cv_for_replicates function
#sampleid: specific sample to plot
#cellid: spcific cell type to plot
cvplot<-function(cvlist,sampleid,cellid){
id1<-which(names(cvlist)==sampleid)
id2<-which(names(cvlist[[id1]])==cellid)
plotdata<-cvlist[[id1]][[id2]]
if(plotdata=='NA'){print(paste0('No ',cellid,' in ', sampleid))}
else{
par(mar = c(8, 8, 8, 8))
boxplot(plotdata,las=2,xlab='Number of cells in replicate',ylab='CV',main=paste0('Sample: ',sampleid,' , ','Cell type: ',cellid))
abline(h=0.1,col='red')
abline(h=0.2,col='blue')}
}

#select sample with the largest number of cell
get_sample<-function(meta){
sample_cell=as.data.frame(rbind(table(meta$sample,meta$celltype)))#sample-celltype-table
sid<-apply(sample_cell,2,function(x){rownames(sample_cell)[which.max(x)]})
snum<-apply(sample_cell,2,function(x){x[which.max(x)]})
return(data.frame(cbind(sid,snum)))
}


get_perm<-function(meta,counts,cvlist,k){
stable<-table(meta$sample,meta$celltype)
pval<-list()
testcv<-list()
for(j in 1:ncol(stable)){
sampleid<-rownames(stable)[which.max(stable[,j])]
cellid<-colnames(stable)[[j]]
num<-names(cvlist[[sampleid]][[cellid]])[[length(names(cvlist[[sampleid]][[cellid]]))]]
num<-as.numeric(num)

rawcounts<-counts[,meta$sample==sampleid&meta$celltype==cellid]  
cvtmp<-list()
for(i in 1:k){
df<-build_rep(rawcounts,num,3)

cvtmp[[i]]<-apply(df,1,get_cv)
}
cvtmp2<-Reduce(cbind,cvtmp)

mincv<-cvlist[[sampleid]][[cellid]][[as.character(num)]]
testcv[[j]]<-cvtmp2[match(names(mincv),rownames(cvtmp2)),]}
names(testcv)<-colnames(stable)
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
