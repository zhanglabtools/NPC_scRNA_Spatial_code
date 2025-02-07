library(tidyverse)

for (pi in c("P1","P3","P5","P6","P7","P9","R3","R6","R7","R8","R9","R10","R12","R13","R14")) {
  usage.file=dir("./important_results/",pattern = paste(pi,"_cNMF.usages",sep = ""))
  usage.df=read.table(paste("./important_results/",usage.file,sep = ""),header = T,sep = "\t",stringsAsFactors = F)
  usage.df=usage.df[,-1]
  colnames(usage.df)=paste(pi,1:dim(usage.df)[2],sep = "_")
  #normalize
  usage.df=usage.df / rowSums(usage.df)
  write.table(usage.df,file = paste(pi,"program.usage.norm.txt",sep = "_"),quote = F,sep = "\t",row.names = T,col.names = T)
  #QC1
  tmpdf1=gather(usage.df,"program","ratio")
  tmpdf1%>%ggplot(aes(x=program,y=ratio))+geom_boxplot(outlier.shape = NA)+geom_jitter(color="red",alpha=0.4)
  ggsave(paste(pi,"program.usage.norm.QC.png",sep = "_"),device = "png",width = 30,height = 16,units = c("cm"))
  #score
  score.file=dir("./important_results/",pattern = paste(pi,"_cNMF.gene_spectra_score",sep = ""))
  score.df=read.table(paste("./important_results/",score.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  score.df=as.data.frame(t(score.df))
  colnames(score.df)=paste(pi,1:dim(score.df)[2],sep = "_")
  
  top30.df=as.data.frame(matrix(nrow = 30,ncol = ncol(score.df)))
  colnames(top30.df)=colnames(score.df)
  for (k in colnames(score.df)) {
    tmpv=score.df[,k]
    names(tmpv)=rownames(score.df)
    top30.df[,k]=names(rev(tail(sort(tmpv),30)))
  }
  write.table(top30.df,file = paste(pi,"program.Zscore.top30gene.txt",sep = "_"),quote = F,sep = "\t",row.names = F,col.names = T)
  score.df$gene=rownames(score.df)
  write.table(score.df,file = paste(pi,"program.Zscore.txt",sep = "_"),quote = F,sep = "\t",row.names = F,col.names = T)
}
############################################################
check.usage=data.frame()
for (pi in c("P1","P3","P5","P6","P7","P9","R3","R6","R7","R8","R9","R10","R12","R13","R14")) {
  usage.file=paste(pi,"_program.usage.norm.txt",sep = "")
  usage.df=read.table(usage.file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  check.usage=rbind(check.usage,as.data.frame(colMeans(usage.df)))
}
colnames(check.usage)=c("mean_ratio")

check.usage$patient_programs=factor(rownames(check.usage),levels = rownames(check.usage))
check.usage%>%ggplot(aes(x=patient_programs,y=mean_ratio))+geom_point()+
  geom_hline(yintercept = 0.03,color="red")+ 
  theme(
    axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)
  )
ggsave("check.usage.1.png",width = 30,height = 16,device = "png",units = "cm")

check.usage=check.usage%>%arrange(mean_ratio)
check.usage$patient_programs=as.character(check.usage$patient_programs)
check.usage$patient_programs=factor(check.usage$patient_programs,levels = check.usage$patient_programs)
check.usage%>%ggplot(aes(x=patient_programs,y=mean_ratio))+geom_point()+
  geom_hline(yintercept = 0.05,color="red")+
  geom_vline(xintercept = 3.5,color="red")+ 
  theme(
    axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)
  )
ggsave("check.usage.2.png",width = 30,height = 16,device = "png",units = "cm") 

maybe.bg=as.character(check.usage$patient_programs[1:3])
