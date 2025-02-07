library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(wesanderson)

all.score.df=data.frame()
all.score.top30.df=data.frame()

for (pi in c("P1","P3","P5","P6","P7","P9","R3","R6","R7","R8","R9","R10","R12","R13","R14")) {
  score.file=paste("./top30_gene/after_qc/",pi,"_program.Zscore.txt",sep = "")
  score.df=read.table(score.file,header = T,sep = "\t",stringsAsFactors = F)
  if (pi=="P1") {all.score.df=score.df}
  if (pi!="P1") {
    all.score.df=all.score.df%>%inner_join(score.df,by="gene")
  }
  
  score.top30.file=paste("./top30_gene/after_qc/",pi,"_program.Zscore.top30gene.txt",sep = "")
  score.top30.df=read.table(score.top30.file,header = T,sep = "\t",stringsAsFactors = F)
  if (pi=="P1") {all.score.top30.df=score.top30.df}
  if (pi!="P1") {
    all.score.top30.df=cbind(all.score.top30.df,score.top30.df)
  }
}

noise_program = c('P3_1','R6_3','R9_5','R8_4','R12_3','R13_5','P1_1','P1_4','R8_1','P9_2','P5_3','R7_3','P5_1','R13_9','R9_1','R14_1','P3_4','R10_3','P1_3',"R13_4")
rownames(all.score.df)=all.score.df$gene
all.score.df$gene=NULL

all.score.df=all.score.df[rowSums(is.na(all.score.df)) == 0,]
all.score.rm.df=all.score.df[,setdiff(colnames(all.score.df),noise_program)] 
all.score.rm.df.cor=cor(all.score.rm.df,method = "pearson")

all.score.rm.df.cor[all.score.rm.df.cor < 0]=0 

colanno=as.data.frame(colnames(all.score.rm.df.cor))
colnames(colanno)="colnames"
colanno$patient=str_replace(colanno$colnames,"_.*","")
rownames(colanno)=colanno$colnames
colanno$colnames=NULL

rowanno=as.data.frame(rownames(all.score.rm.df.cor))
colnames(rowanno)="rownames"
rowanno$patient=str_replace(rowanno$rownames,"_.*","")
rownames(rowanno)=rowanno$rownames
rowanno$rownames=NULL
#write.csv(rowanno,'rowanno.csv')

rowanno1 = read.csv('rowanno.csv', row.names = 'X')

#color_v=brewer.pal(2, "Set3")
#names(color_v)=c("P1","P3","P5","P6","P7","P9","R3","R6","R7","R8","R9","R10","R12","R13","R14")
color_v = c('#5f82b5','#D76364','#BEBADA')
names(color_v)=c('P','R')
ann_colors = list(patient = color_v)

#pushViewport(viewport(gp = gpar(fontfamily = "Arial")))
p = pheatmap(all.score.rm.df.cor,cluster_rows = T,cluster_cols = T,
         clustering_method = "ward.D2", 
         show_colnames = F,
         treeheight_row=0,treeheight_col=0,
         border_color=NA,
         annotation_row = rowanno1, annotation_col = rowanno1,
         annotation_names_row = F,annotation_names_col = F,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#f9fbfd","#eae873","#ecaa37", "#b3403b"))(50),
         fontsize_row=8,
         fontsize = 15,
         width = 11.5,height = 9,
         filename = 'program_pearson_cor.ward.D2.heatmap.pdf'
)
p
all.score.top30.rm.df=all.score.top30.df[,setdiff(colnames(all.score.top30.df),c("P5_3","P5_4","P5_5","P1_3","P1_6","P1_2","P1_7"))]#在质控这一步检测出来的噪声
write.table(all.score.top30.rm.df,file = "program_top30gene.txt",quote = F,sep = "\t",row.names = F,col.names = T)

#实际分析中的一些经???
#相关性聚类这一步需要多做几次。根据每个program的大致功能，如果发现另一种功能的program聚到某一种meta模块里面，这时可以将乱入的program删掉，再做一次相关性聚??????
#如果我们认定应该属于同一个meta模块的program分散在两个地方，可以试试调整聚类方法(参数clustering_method)，或者像"_1""_2"这样定义
