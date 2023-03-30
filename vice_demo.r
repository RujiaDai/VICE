#Demo for VICE package
load('cv_sc.r')
meta<-read.csv('metadata.csv',header=T,sep=',') #load metadata for cells, must include information: sample id, cell identity
counts<-read.csv('matrix.csv',header=T,sep=',',row.names=1)# load count matrix, genes by cells
counts<-counts[match(meta$sample_name,rownames(counts)),]# make sure count matrix and metadata are in same order



meta_class<-meta[,c('external_donor_name_id','class_label')] # extract sample id and cell identity from metadata
colnames(meta_class)<-c('sample','celltype')


cvlist1<-get_cv_for_replicates(meta_class,counts)#calculate CV value for each gene, each cell type, in each sample

stable1<-get_sample(meta_class)#get the sample id which has the largest cell numbers


perm1<-get_perm(meta_class,counts,cvlist1,100)# 100 times permutation test for the CV calculated in the sample in stable1


save(oerm1,cvlist1,stable1,file='aba_mouse_class_cvlist.RData')
