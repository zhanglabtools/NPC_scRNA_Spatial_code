library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(xlsx)

hsets <- read.gmt("hallmark_cancersea.gmt")
enrich.result=data.frame()

program_top30=read.table("program_top30gene.txt",header = T,sep = "\t",stringsAsFactors = F)
for (i in 1:dim(program_top30)[2]) {
  tmp <- enricher(program_top30[,i], TERM2GENE = hsets)
  if (is.null(tmp)) {
    next
  }
  tmp_result <- tmp@result
  tmp1=head(tmp_result)
  tmp1$program=colnames(program_top30)[i]
  rownames(tmp1)=NULL
  enrich.result=rbind(enrich.result,tmp1)
}
write.csv(enrich.result,file = "program_top30gene_anno.csv",row.names = F)
