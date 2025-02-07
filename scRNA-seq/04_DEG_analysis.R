library(Matrix)
library(data.table)
library(Seurat)
#library(SeuratDisk)
library(dplyr)
#library(gt)
library(EnhancedVolcano)
#library(igraph)
library(RColorBrewer)


data = Seurat::Read10X('E:/CNV annot malignant/Malignant', gene.column = 1)
Data = CreateSeuratObject(data, min.cells = 3, min.genes = 200)
write.csv(rownames(Data@assays$RNA@counts),'total_genes.csv')
Meta_data = read.table('E:/CNV annot malignant/Mali_meta.csv', header=TRUE, sep=',', row.names='X')
#Data$patients = Meta_data$V1
Data$type = Meta_data$type
Data$patient = Meta_data$PID

Data <- Seurat::NormalizeData(Data, verbose = FALSE) %>%
  Seurat::FindVariableFeatures(., nfeatures = 3000) %>%
  Seurat::ScaleData(.,) %>%
  Seurat::RunPCA(., verbose = FALSE) %>%
  Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)

Seurat::Idents(object = Data) <- Data$type
cluster_markers_all <- Seurat::FindAllMarkers(object = Data, logfc.threshold = -Inf,
                                              slot = "data",
                                              min.diff.pct = 0,
                                              verbose = TRUE, 
                                              only.pos = TRUE)
write.csv(cluster_markers_all,'DEGs_Malignant_GSEA.csv')
#topgene <- cluster_markers_all$gene[cluster_markers_all$p_val_adj<1e-100]
cluster_markers_all[cluster_markers_all$cluster=='P', ]$avg_log2FC = cluster_markers_all[cluster_markers_all$cluster=='P', ]$avg_log2FC*(-1)
###new version###
library(ImageGP)
library(ggplot2)
library(ggpubr)
library(egg)
library(ggrepel)
PT_index = Meta_data$type=='P'
RT_index = Meta_data$type=='R'
Gene_expr_mean_PT = rowMeans(Data@assays[["RNA"]][,PT_index])[cluster_markers_all$gene]
Gene_expr_mean_RT = rowMeans(Data@assays[["RNA"]][,RT_index])[cluster_markers_all$gene]
cluster_markers_all['PT'] = Gene_expr_mean_PT
cluster_markers_all['RT'] = Gene_expr_mean_RT
cluster_markers_all$level <- ifelse(cluster_markers_all$p_val_adj<0.05, 
                         ifelse(cluster_markers_all$avg_log2FC>=0.5, "Up", 
                                ifelse(cluster_markers_all$avg_log2FC<=-0.5, "Down", "NoSig")),"NoSig")
head(cluster_markers_all)
#定义横纵坐标变量，用的是前面计算的样本平均值
                    size_variable=1.2,
                    color_variable = "level",
                    title ="PT vs RT", 
                    color_variable_order = c("NoSig","Up", "Down"),
                    manual_color_vector = c("grey","#d76364","#5f82b5")) +
                    labs(x = "RT", y = "PT")
p
ggsave('volcano_plot_45.png', dpi=600, width = 10, height = 8)

cluster_markers_all['gene'] = rownames(cluster_markers_all)
cluster_markers_all$label =""
cluster_markers_all <- cluster_markers_all[order(cluster_markers_all$p_val_adj),]
up.genes <- head(cluster_markers_all$gene[which(cluster_markers_all$level=="Up")],10)
down.genes <- head(cluster_markers_all$gene[which(cluster_markers_all$level=="Down")],10)
top10genes <- c(as.character(up.genes), as.character(down.genes))
cluster_markers_all$label[match(top10genes,cluster_markers_all$gene)] <- top10genes

p1 = p +  geom_text_repel(data=cluster_markers_all, aes(label= label), color="black", size=5, fontface="italic",
                     arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                     point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, 
                     max.iter = 3e3)
p1
