library(Seurat)
library(monocle3)
library(ggpubr)
library(Matrix)
library(ComplexHeatmap)
library(igraph)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(patchwork)
library(colorRamps)
library(cowplot)
library(VGAM)
rm(list=ls())
empty_theme <- theme_void() + theme(legend.position="none")

expr_matrix <- read.csv('E:\\CAF\\save\\CAF_data.csv',header=TRUE,row.names = 'X')
expr_matrix <- t(expr_matrix)

HVGs = read.csv('E:\\CAF\\save\\CAF_hvgs3000.csv')
Meta_data = read.table('E:\\CAF\\save\\CAF_meta.csv',sep=',',header=TRUE,row.names = 'X')
#Meta_data = Meta_data[colnames(expr_matrix), ]
data_seurat = CreateSeuratObject(counts = expr_matrix, project = "CAF", min.cells = 3, min.features = 1, meta.data=Meta_data)
data_seurat <- NormalizeData(data_seurat, normalization.method = "LogNormalize", scale.factor = 1000)
#data_seurat <- AddMetaData(data_seurat, metadata = PID)

#data <- GetAssayData(data_seurat, assay = 'RNA', slot = 'counts')
data = expr_matrix
cell_metadata <- data_seurat@meta.data
allgenes <- row.names(data_seurat@assays$RNA@data)
gene_annotation <- data.frame(gene=allgenes,gene_short_name=allgenes, use_for_ordering=F)
row.names(gene_annotation) <- allgenes
data = data[allgenes, ]

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50, norm_method="none", use_genes = intersect(HVGs$X0,rownames(rowData(cds))))
harmony_embed = read.csv('E:\\CAF\\save\\X_pca_harmony.csv',row.names = 'X')
harmony_umap = read.csv('E:\\CAF\\save\\X_umap.csv',row.names = 'X')
cds@int_colData$reducedDims$UMAP <- as.matrix(harmony_umap)
cds@int_colData$reducedDims$PCA <- as.matrix(harmony_embed)


cds@clusters@listData$UMAP$partitions <- Meta_data$celltype
names(cds@clusters@listData$UMAP$partitions) <- rownames(Meta_data)
cds@clusters@listData$UMAP$clusters <- Meta_data$celltype
names(cds@clusters@listData$UMAP$clusters) <- rownames(Meta_data)
#colnames(cds) = rownames(Meta_data)
cds <- learn_graph(cds, use_partition=FALSE, close_loop=FALSE)
cds <- order_cells(cds, reduction_method="UMAP")
CAF_pseudotime = cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
write.csv(CAF_pseudotime, 'CAF_pseudotime.csv')
my_color = c('#365180','#ba86b4','#f4a675','#fcdb72')
names(my_color) = c('CAF MCAM','CAF GREM1','CAF CFD','CAF SFRP4')
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="type")
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",label_groups_by_cluster = FALSE, label_cell_groups = FALSE, trajectory_graph_color = "#0a0a0a",
                 label_leaves = FALSE, cell_size = 0.8, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,trajectory_graph_segment_size = 1.25,)
p4 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype",label_groups_by_cluster = FALSE, label_cell_groups = FALSE, trajectory_graph_color = "#0a0a0a",
                 label_leaves = FALSE, cell_size = 0.8, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,trajectory_graph_segment_size = 1.25,)+ 
  scale_color_manual(values = my_color) 
p3 + p4

Track_genes <- graph_test(cds, neighbor_graph = 'principal_graph', cores=1)
genes <- row.names(subset(Track_genes, q_value < 0.05 & morans_I > 0.2))
match_genes <- match(genes,rownames(rowData(cds)))
match_genes <- na.omit(match_genes)
pt.matrix <- exprs(cds)[match_genes, order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- rownames(rowData(cds))[match_genes]

tiff('heatmap.tiff', res=600, width = 2580, height = 2000, compression='lzw')
kclust <- kmeans(pt.matrix, 3)
split <- paste0("Cluster\n", kclust$cluster)
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=6),rev(brewer.pal(6, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 2.5),
  km = 3,
  #heatmap_height = nrow(pt.matrix)*unit(0.15, 'cm'),
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,)
#split=split)
#font_family = "arial"
#htkm@column_names_param$gp$fontfamily = font_family  # ??????
#htkm@row_names_param$gp$fontfamily = font_family   # ??????
#htkm@row_title_param$gp$fontfamily = font_family  # ?б???p@column_title_param$gp$fontfamily = font_family  # ?б??⣨????????
htkm
dev.off()
