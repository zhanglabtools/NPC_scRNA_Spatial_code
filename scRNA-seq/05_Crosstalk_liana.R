library(tidyverse)
library(magrittr)
library(ggplot2)
library(liana)
library(Seurat)
show_resources()
show_methods()

expression <- Seurat::Read10X('E:/NPC_data/total_data', gene.column = 1)
sample_info <- read.table('E:/NPC_data/Mali_meta.csv', header=TRUE, sep=',', row.names='X')
# create seurat object
seurat_object <- Seurat::CreateAssayObject(counts = expm1(t(expression))) %>%
  Seurat::CreateSeuratObject(., meta.data = sample_info) %>%
  Seurat::NormalizeData()

# set cell identity to cell type
Idents(seurat_object) <- seurat_object@meta.data$cell_type
liana_results <- liana_wrap(seurat_object) %>%
  liana_aggregate()

liana_trunc <- liana_test %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.05) # note that these pvals are already corrected

write.csv(liana_results, 'liana_results.csv')
write.csv(liana_trunc, 'liana_trunc.csv')

