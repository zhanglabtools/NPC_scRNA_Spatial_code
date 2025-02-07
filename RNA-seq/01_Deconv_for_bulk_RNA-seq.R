library(Seurat)
suppressWarnings(library(BayesPrism))

expr_matrix <- read.csv('D:\\DataStore\\alldata.csv', row.names = 'X')
sc_annotation = read.csv('D:\\deconv_NPC\\alldata_meta.csv', row.names = 'X')
sc_counts <- as.matrix(expr_matrix)
bulk_matrix <- read.table('D:\\deconv_NPC\\NPC_不良预后\\GSE102349_NPC_mRNA_processed.txt', sep='\t', header = TRUE, row.names = 1)
bulk_matrix <- t(as.matrix(bulk_matrix))
sort(table(sc_annotation$celltype1))
plot.cor.phi(input=sc_counts,
             input.labels=sc_annotation$celltype1,
             title="cell type correlation",
             #specify pdf.prefix if need to output to pdf
             #pdf.prefix="gbm.cor.ct",
             cexRow=0.5, cexCol=0.5,
)

sc.stat <- plot.scRNA.outlier(
  input=sc_counts, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=sc_annotation$celltype1,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting.
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)
bk.stat <- plot.bulk.outlier(
  bulk.input=bulk_matrix,#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=sc_counts, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=sc_annotation$celltype1,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)
sc.dat.filtered <- cleanup.genes(input=sc_counts,
                                 input.type="count.matrix",
                                 species="hs",
                                 gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"),
                                 exp.cells=5)
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bulk_matrix
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,gene.type = "protein_coding")
diff.exp.stat <- get.exp.stat(sc.dat=sc_counts[,colSums(sc_counts>0)>3],# filter genes to reduce memory use
                              cell.type.labels=sc_annotation$celltype1,
                              cell.state.labels=paste(sc_annotation$PID,sc_annotation$celltype1),
                              psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value n.cores=1 #number of threads
)
sc.dat.filtered.pc.sig <- select.marker(sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,                                                           
                                         lfc.min=0.1)
myPrism <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=bulk_matrix,
  input.type="count.matrix",
  cell.type.labels = sc_annotation$celltype1,
  cell.state.labels = sc_annotation$sub_class,
  key="Malignant",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
bp.res <- run.prism(prism = myPrism, n.cores=20)
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
theta1 <- get.fraction (bp=bp.res,
                       which.theta="first",
                       state.or.type="state")
write.csv(theta, 'theta.csv')
write.csv(theta1, 'theta1.csv')
save.image('celltype.rdata')
