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
expr_matrix <- Seurat::Read10X('E:\\CD8 T\\monocle3\\CD8T_new\\CD8T_cell',gene.column = 1)

HVGs = read.csv('E:\\CD8 T\\monocle3\\CD8T_new\\CD8T_hvgs3000.csv')
Meta_data1 = read.table('E:\\CD8 T\\monocle3\\CD8T_new\\CD8T_meta.csv',sep=',',header=TRUE,row.names = 'X')
Meta_data = read.table('E:\\CD8 T\\monocle3\\CD8T_new\\CD8T_meta2.csv',sep=',',header=TRUE,row.names = 'X')
Meta_data$sub_class1 = Meta_data1$sub_class1
data_seurat = CreateSeuratObject(counts = expr_matrix, project = "CD8T", min.cells = 3, min.features = 200, meta.data=Meta_data)
data_seurat <- NormalizeData(data_seurat, normalization.method = "LogNormalize", scale.factor = 1000)
#data_seurat <- AddMetaData(data_seurat, metadata = PID)

#data <- GetAssayData(data_seurat, assay = 'RNA', slot = 'counts')
data = data_seurat@assays$RNA@data
data = matrix(data, nrow=dim(data)[1], ncol=dim(data)[2])
cell_metadata <- data_seurat@meta.data
allgenes <- row.names(data_seurat@assays$RNA@data)
gene_annotation <- data.frame(gene=allgenes,gene_short_name=allgenes, use_for_ordering=F)
row.names(gene_annotation) <- allgenes

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50, norm_method="none", use_genes = HVGs$X0)
harmony_embed = read.csv('E:\\CD8 T\\monocle3\\CD8T_new\\CD8T_harmony_pca.csv',row.names = 'X')
harmony_umap = read.csv('E:\\CD8 T\\monocle3\\CD8T_new\\CD8T_umap.csv',row.names = 'X')
cds@int_colData$reducedDims$UMAP <- as.matrix(harmony_umap)
cds@int_colData$reducedDims$PCA <- as.matrix(harmony_embed)


cds@clusters@listData$UMAP$partitions <- Meta_data$sub_class1
names(cds@clusters@listData$UMAP$partitions) <- rownames(Meta_data)
cds@clusters@listData$UMAP$clusters <- Meta_data$sub_class1
names(cds@clusters@listData$UMAP$clusters) <- rownames(Meta_data)
#colnames(cds) = rownames(Meta_data)
cds <- learn_graph(cds, use_partition=FALSE, close_loop=FALSE)
cds <- order_cells(cds, reduction_method="UMAP")
CD8T_pseudotime = cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
write.csv(CD8T_pseudotime, 'CD8T_pseudotime.csv')

cell_proj_tree = cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_tree"]]

#degrees <- degree(cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_tree"]])
#write.csv(degrees,'./graph_save/degrees.csv')
#closest_vertex <- cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_closest_vertex"]]
#closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#write.csv(closest_vertex,'./graph_save/closest_vertex.csv')
#vertices_means = read.csv('E:/zhongshan/code_04_17/pseudotime trajectory/monocle3/CD8T_new/graph_save/vertices_means.csv',row.names='X')

reduction_method = 'UMAP'
ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>% 
as.data.frame() %>% 
dplyr::select(prin_graph_dim_1 = X0, prin_graph_dim_2 = X1) %>% 
dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

dp_mst <- cds@principal_graph[[reduction_method]]
edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select(source = "from", 
target = "to") %>% dplyr::left_join(ica_space_df %>% 
dplyr::select(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", 
source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% 
dplyr::left_join(ica_space_df %>% 
dplyr::select(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

plot_a <- ggplot(ica_space_df, aes(x=prin_graph_dim_1, y=prin_graph_dim_2)) +
  geom_point(color='black', size=0.2) + 
  xlim(0, 12.0) +
  ylim(-3, 10.0) + 
  geom_label(data=ica_space_df,
             aes(prin_graph_dim_1,prin_graph_dim_2,label=sample_name, colour='blue'))+
  theme(legend.position="none")
p1 <- plot_cells(cds, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, label_leaves = FALSE, cell_size = 0.2,alpha=0.8,
                 label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,trajectory_graph_segment_size = 1.25,trajectory_graph_color = "#0a0a0a",)+
  scale_color_manual(values=c('#b2df8a', '#fdbf6f', '#a6cee3'))
p1
plot_grid(plot_a, p1 , align = "h", labels = c('Centroids (graph nodes)', 'Cells'))


p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="type")
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",label_groups_by_cluster = FALSE, label_cell_groups = FALSE, trajectory_graph_color = "#0a0a0a",
                 label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,trajectory_graph_segment_size = 1,)
p4 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="sub_class",label_groups_by_cluster = FALSE, label_cell_groups = FALSE, trajectory_graph_color = "#0a0a0a",
                 label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,trajectory_graph_segment_size = 1,)+
                  scale_color_manual(values=c('#f8766d','#cd9600', '#7cae00','#00be67','#00bfc4', '#ffa257','#c77bff','#ff59ca'))
p1+p3
p3+p4

umap_branch = read.csv('E:\\CD8 T\\monocle3\\CD8T_new\\graph_save\\Umap_branch.csv', row.names = 'X')
colData(cds)$cell_fate = umap_branch$cellfate
cds_cyto <- cds[, colData(cds)$cell_fate %in% c(1,3)]
cds_exha <- cds[, colData(cds)$cell_fate %in% c(1,2)]
p5 <-plot_cells(cds_cyto, reduction_method="UMAP", color_cells_by="pseudotime",label_groups_by_cluster = FALSE, label_cell_groups = FALSE, trajectory_graph_color = "#0a0a0a",
           label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,trajectory_graph_segment_size = 1.25,)
p6 <- plot_cells(cds_exha, reduction_method="UMAP", color_cells_by="pseudotime",label_groups_by_cluster = FALSE, label_cell_groups = FALSE, trajectory_graph_color = "#0a0a0a",
           label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,trajectory_graph_segment_size = 1.25,)
p5+p6
Track_genes1 <- graph_test(cds_cyto, neighbor_graph = 'principal_graph', cores=12)###windows??֧?ֲ???
#write.csv(Track_genes1, 'Track_genes_cyto.csv')
Track_genes2 <- graph_test(cds_exha, neighbor_graph = 'principal_graph', cores=12)###windows??֧?ֲ???
#write.csv(Track_genes2, 'Track_genes_exha.csv')


Track_genes <- graph_test(cds, neighbor_graph = 'principal_graph', cores=12)###windows??֧?ֲ???
#write.csv(Track_genes, 'Track_genes.csv')
#Track_genes <- read.csv('Track_genes.csv', row.name='X')
deg_ids <- row.names(subset(Track_genes, q_value==0))

###Gene module
gene_module_df <- find_gene_modules(cds[deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=partitions(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)

module_dendro <- hclust(dist(agg_mat))
#gene_module_df$module <- factor(gene_module_df$module, levels = row.names(agg_mat)[module_dendro$order])
write.csv(gene_module_df, 'gene_module_df.csv')
plot_cells(cds,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE)

select_genes = c('CCR7','TCF7','GZMA','GZMH','GZMK','GNLY','HAVCR2','TIGIT','CTLA4')
select_lineage_cds <- cds[rowData(cds)$gene_short_name %in% select_genes,
                       ]
monocle3::plot_genes_in_pseudotime(select_lineage_cds,
                         color_cells_by="sub_class",
                         min_expr=0.05, ncol = 2)

Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#????????????ͼ
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="sub_class", 
                         min_expr=0.5, ncol = 2)

####??P_R?????Ļ???
P_cds <- cds[, colData(cds)$type %in% c("P")]
R_cds <- cds[, colData(cds)$type %in% c("R")]
gene_list = c('CCR7','TCF7','GZMA','GZMH','GZMK','GNLY','HAVCR2','TIGIT','CTLA4')
gene_list = c('BCL11B','FOXP1','KLF10','HOPX','NR4A1','TOX')
gene_list = c('BTLA','CTLA4','HAVCR2','LAG3','PDCD1','TIGIT')
gene_list = c('CD8A','CD8B')
p10<-plot_genes_in_pseudotime(P_cds[gene_list,], color_cells_by="sub_class", min_expr=0.05, ncol = 2)
p11<-plot_genes_in_pseudotime(R_cds[gene_list,], color_cells_by="sub_class", min_expr=0.05, ncol = 2)
p10+p11
p4 <- plot_cells(P_cds, reduction_method="UMAP", color_cells_by="sub_class",label_groups_by_cluster = FALSE, label_cell_groups = FALSE, 
                 label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,
                 trajectory_graph_segment_size = 1.25,trajectory_graph_color = "#0a0a0a")+
                 scale_color_manual(values=c('#f8766d','#cd9600', '#7cae00','#00be67','#00bfc4', '#ffa257','#c77bff','#ff59ca'))
p5 <- plot_cells(R_cds, reduction_method="UMAP", color_cells_by="sub_class", label_cell_groups=FALSE,label_groups_by_cluster=FALSE, 
                 label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,
                 trajectory_graph_segment_size = 1.25,trajectory_graph_color = "#0a0a0a")+
                 scale_color_manual(values=c('#f8766d','#cd9600', '#7cae00','#00be67','#00bfc4', '#ffa257','#c77bff','#ff59ca'))
p4+p5

p6 <- plot_cells(P_cds, reduction_method="UMAP", genes=c("CCR7"),label_groups_by_cluster = FALSE, label_cell_groups = FALSE, 
                 label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,
                 trajectory_graph_segment_size = 1.25,trajectory_graph_color = "#0a0a0a")
p7 <- plot_cells(R_cds, reduction_method="UMAP", genes=c("GNLY"), label_cell_groups=FALSE,label_groups_by_cluster=FALSE, 
                 label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,
                 trajectory_graph_segment_size = 1.25,trajectory_graph_color = "#0a0a0a")
p8 <- plot_cells(R_cds, reduction_method="UMAP", genes=c("TIGIT"), label_cell_groups=FALSE,label_groups_by_cluster=FALSE, 
                 label_leaves = FALSE, cell_size = 0.2, label_branch_points = FALSE, label_roots=FALSE, group_label_size=5,
                 trajectory_graph_segment_size = 1.25,trajectory_graph_color = "#0a0a0a")
p6+p7

###cell density
df <- as.data.frame(cds@colData@listData)
df$pseudotime <- CD8T_pseudotime
p1 = ggplot(df[df$type %in% c('P'), ], aes(pseudotime, colour=orig.ident, fill=orig.ident)) + geom_density(bw=0.5,size=1,alpha=0.5,color='#5f82b5',fill='#5f82b5')+ylim(0,0.2)+theme_classic2()
p2 = ggplot(df[df$type %in% c('R'), ], aes(pseudotime, colour=orig.ident, fill=orig.ident)) + geom_density(bw=0.5,size=1,alpha=0.5,color='#d76364',fill='#d76364')+ylim(0,0.2)+theme_classic2()
p1+p2

Track_genes <- read.csv('Track_genes.csv', row.name='X')
genes <- row.names(subset(Track_genes, q_value < 0.01 & morans_I > 0))
cds_subset <- cds[rowData(cds)$gene_short_name %in% genes,]
gene_fits_pseudotime <- fit_models(cds_subset, model_formula_str = "~pseudotime")
fit_coefs_pseudotime <- coefficient_table(gene_fits_pseudotime)
fit_coefs_pseudotime_terms <- fit_coefs_pseudotime %>% filter(term == "pseudotime")
fit_coefs_pseudotime_terms = fit_coefs_pseudotime_terms %>% filter (q_value < 0.05) %>%select(gene, term,p_value, q_value, estimate)
df = cbind(fit_coefs_pseudotime_terms$gene,fit_coefs_pseudotime_terms$p_value, fit_coefs_pseudotime_terms$q_value, fit_coefs_pseudotime_terms$term)
colnames(df) = c('gene','p value', 'q value', 'term')
write.csv(df,'gene_fits_pseudotime.csv')
gene_fits_cluster <- fit_models(cds_subset, model_formula_str = "~sub_class1")
fit_coefs_cluster <- coefficient_table(gene_fits_cluster)
fit_coefs_cluster_terms <- fit_coefs_cluster %>% filter(term != "(Intercept)")
fit_coefs_cluster_terms = fit_coefs_cluster_terms %>% filter (q_value < 0.05) %>%select(gene, term,p_value, q_value, estimate)
df = cbind(fit_coefs_cluster_terms$gene,fit_coefs_cluster_terms$p_value, fit_coefs_cluster_terms$q_value, fit_coefs_cluster_terms$term)
colnames(df) = c('gene','p value', 'q value', 'term')
write.csv(df,'gene_fits_cluster.csv')

####heatmap
Track_genes <- read.csv('Track_genes.csv', row.name='X')
genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.1))
#genes_scor <- read.csv('E:\\中大\\22_04_17\\Pseudotime trajectory\\CD8 T\\monocle3\\R code\\gene_sel_scorpius.csv')
#genes <- genes_scor$x
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
  row_names_gp                 = gpar(fontsize = 5.5),
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

#####Plot branches
Track_genes1 = read.csv('Track_genes_cyto.csv', row.names = 'X')
Track_genes2 = read.csv('Track_genes_exha.csv', row.names = 'X')
genes1 <- row.names(subset(Track_genes1, q_value < 0.05 & morans_I > 0.1))
genes2 <- row.names(subset(Track_genes2, q_value < 0.05 & morans_I > 0.1))
gene = intersect(genes1, genes2)
match_genes1 <- match(genes,rownames(rowData(cds_cyto)))
match_genes1 <- na.omit(match_genes1)
pt.matrix1 <- exprs(cds_cyto)[match_genes1, order(pseudotime(cds_cyto),decreasing =TRUE)]
pt.matrix1 <- t(apply(pt.matrix1,1,function(x){smooth.spline(x,df=3)$y}))
#pt.matrix1 <- t(apply(pt.matrix1,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix1) <- rownames(rowData(cds_cyto))[match_genes1]

match_genes2 <- match(genes,rownames(rowData(cds_exha)))
match_genes2 <- na.omit(match_genes2)
pt.matrix2 <- exprs(cds_exha)[match_genes2, order(pseudotime(cds_exha))]
pt.matrix2 <- t(apply(pt.matrix2,1,function(x){smooth.spline(x,df=3)$y}))
#pt.matrix2 <- t(apply(pt.matrix2,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix2) <- rownames(rowData(cds_exha))[match_genes2]

pt.matrix_all = cbind(pt.matrix1,pt.matrix2)
kclust <- kmeans(pt.matrix_all, 5)
split <- paste0("Cluster\n", kclust$cluster)
tiff('heatmap_cyto.tiff', res=600, width = 2580, height = 2000, compression='lzw')
htkm <- Heatmap(
  pt.matrix1,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-v,to=2,length=6),rev(brewer.pal(6, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 2),
  #km = 5,
  #heatmap_height = nrow(pt.matrix)*unit(0.15, 'cm'),
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
split=split)
#font_family = "arial"
#htkm@column_names_param$gp$fontfamily = font_family  # ??????
#htkm@row_names_param$gp$fontfamily = font_family   # ??????
#htkm@row_title_param$gp$fontfamily = font_family  # ?б???p@column_title_param$gp$fontfamily = font_family  # ?б??⣨????????
htkm
dev.off()

tiff('heatmap_exha.tiff', res=600, width = 2580, height = 2000, compression='lzw')
#kclust <- kmeans(pt.matrix, 3)
split <- paste0("Cluster\n", kclust$cluster)
htkm <- Heatmap(
  pt.matrix2,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=6),rev(brewer.pal(6, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 2),
  #km = 5,
  #heatmap_height = nrow(pt.matrix)*unit(0.15, 'cm'),
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  split=split)
#font_family = "arial"
#htkm@column_names_param$gp$fontfamily = font_family  # ??????
#htkm@row_names_param$gp$fontfamily = font_family   # ??????
#htkm@row_title_param$gp$fontfamily = font_family  # ?б???p@column_title_param$gp$fontfamily = font_family  # ?б??⣨????????
htkm
dev.off()

scale_max = 2
scale_min = -2
heatmap_matrix = cbind(pt.matrix1,pt.matrix2)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1, sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix), center = TRUE))
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
col_gap_ind = dim(pt.matrix1)[2]
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
exp_rng <- range(heatmap_matrix)
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
hclust_method = c( "ward.D", "ward.D2", "single", "complete", "average","mcquitty","median","centroid")
hclust = hclust_method[1]
save_place=paste(hclust, '.tiff',sep='')
tiff(save_place, res=600, width = 2580, height = 2000, compression='lzw')
ph <- pheatmap(heatmap_matrix, cluster_cols = FALSE, use_raster=T,
               cluster_rows = TRUE, show_rownames = T, show_colnames = F, 
               clustering_distance_rows = row_dist, clustering_method = hclust, 
               gaps_col = col_gap_ind,fontsize=1, treeheight_row = 20,
               cutree_rows = 6,  filename = NA, 
               breaks = bks, color = hmcols)
ph
dev.off()
