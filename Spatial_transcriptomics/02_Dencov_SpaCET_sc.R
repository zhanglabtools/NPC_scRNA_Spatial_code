library(Seurat)
library(SpaCET)
#setwd("D:/Users/qlshen/zhongshan/Figure/Malignant/deconv")

expr_matrix <- read.csv('D:\\DataStore\\alldata.csv', row.names = 'X')
sc_annotation = read.csv('D:\\DataStore\\alldata_meta.csv', row.names = 'X')

file_list = list.files('D:\\datasets\\space\\')
sc_counts <- t(as.matrix(expr_matrix))

sc_lineageTree = list("Myeloid","B","Malignant","Plasmacyte","Mast","pDC","Fibroblast","Neutrophil","Endothelial","CD4T","CD8T",'NK',"Treg","Tcycling")
names(sc_lineageTree) = c("Myeloid","B","Malignant","Plasmacyte","Mast","pDC","Fibroblast","Neutrophil","Endothelial","CD4T","CD8T",'NK',"Treg","Tcycling")
sc_lineageTree$B = c('B Memory', 'B IGHM', 'B FCRL4', 'B NEIL1')
sc_lineageTree$Myeloid = c("Macro TNFRSF9","Macro SPP1","DC2 CD1A","DC1","DC2 PLAC8","DC3",
                           "Myeloid cycling","Mono FCN1","Macro CCL5","Macro MSR1","Macro ISG15")
sc_lineageTree$Fibroblast = c("MCAM myCAF","GREM1 myCAF","SFRP4 iCAF","CFD iCAF")
sc_lineageTree$Treg = c("Treg CXCR4","Treg CTLA4")
sc_lineageTree$CD8T = c("CD8 PNISR","CD8 FOSB","CD8 TIGIT","CD8 MZT2B","CD8 HAVCR2","CD8 TCF7","CD8 ZNF683")
sc_lineageTree$CD4T = c("CD4 MAF","CD4 IL7R","CD4 ISG")
sc_lineageTree$NK = c("NK CD16" ,"NK ISG","NK CD81","NK XCL1")


for (i in 1:15){
  patient_name = file_list[i]
  path = sprintf('D:\\Users\\qlshen\\zhongshan\\datasets\\space\\%s', patient_name)
  st_data = Load10X_Spatial(path)
  #counts <- Seurat::Read10X(path, gene.column = 1)
  counts <- st_data@assays[["Spatial"]]@counts
  counts <- as.matrix(counts)
  spotCoordinates <- read.csv(sprintf('D:\\Users\\qlshen\\zhongshan\\datasets\\space\\%s\\spatial\\tissue_positions_list.csv', patient_name),header = FALSE, row.names = 'V1')
  spotCoordinates <- spotCoordinates[, 4:5]
  spotCoordinates <- spotCoordinates[colnames(counts), ]
  colnames(spotCoordinates) = c('X','Y')
  
  SpaCET_obj <- create.SpaCET.object(
    counts=counts,
    spotCoordinates=spotCoordinates,
    imagePath=NA,
    platform = "Visium"
  )
  
  SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
    SpaCET_obj, 
    sc_counts=sc_counts, 
    sc_annotation=sc_annotation, 
    sc_lineageTree=sc_lineageTree, 
    sc_nCellEachLineage = 100,
    coreNo=1
  )
  save_place = sprintf('D:\\deconv\\SpaCET\\save_100\\%s_propMat.csv', patient_name)
  write.csv(SpaCET_obj@results[["deconvolution"]][["propMat"]], save_place)
}
