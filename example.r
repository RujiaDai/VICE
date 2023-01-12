#Example for VICE
#ABA mouse data


setwd('/mnt/compute/Groups/LiuLab/User/Rdai/data/aba_mouse/')
meta<-read.csv('metadata.csv',header=T,sep=',')
counts<-read.csv('matrix.csv',header=T,sep=',',row.names=1)
counts<-counts[match(meta$sample_name,rownames(counts)),]

#counts<-t(counts)

#two resolution levels
meta_class<-meta[,c('external_donor_name_id','class_label')]
colnames(meta_class)<-c('sample','celltype')

meta_subclass<-meta[,c('external_donor_name_id','subclass_label')]
colnames(meta_subclass)<-c('sample','celltype')

cvlist1<-get_cv_for_replicates(meta_class,counts)
cvlist2<-get_cv_for_replicatesmeta_subclass,counts)

#plot
cvplot(cvlist2,'373','L4/5 IT CTX')
cvplot(cvlist2,'373','L5 PPP')

