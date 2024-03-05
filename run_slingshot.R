rm(list = ls())
gc()
library(anndata)
library(SingleCellExperiment)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(slingshot)
library(umap)
library(SeuratDisk)

setwd("/localscratch/ziqi/Prostate_Cancer")
X_castrated <- as(Matrix::readMM("Cell_Ranger_output/X_castrated.mtx"), "dgCMatrix")
meta_castrated <- read.csv("Cell_Ranger_output/meta_castrated.csv", row.names = 1, sep = "\t")
xumap_castrated <- read.csv("Cell_Ranger_output/castrated_umap.csv", row.names = 1, sep = "\t")
gene_castrated <- read.csv("Cell_Ranger_output/gene_castrated.csv", row.names = 1)
# load seurat integrated results
seurat_castrated <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate_castrated.h5Seurat")

rownames(X_castrated) <- rownames(meta_castrated)
colnames(X_castrated) <- rownames(gene_castrated)
X_castrated <- X_castrated[(meta_castrated$annot.transfer == "Basal") | (meta_castrated$annot.transfer == "Luminal") | (meta_castrated$annot.transfer == "Mesenchymal"),]
xumap_castrated <- xumap_castrated[(meta_castrated$annot.transfer == "Basal") | (meta_castrated$annot.transfer == "Luminal") | (meta_castrated$annot.transfer == "Mesenchymal"),]
meta_castrated <- meta_castrated[(meta_castrated$annot.transfer == "Basal") | (meta_castrated$annot.transfer == "Luminal") | (meta_castrated$annot.transfer == "Mesenchymal"),]
seurat_castrated <- seurat_castrated[,rownames(meta_castrated)]
# integrated embedding
xpca_castrated_integrated <- seurat_castrated@reductions$pca@cell.embeddings
xumap_castrated_integrated <- seurat_castrated@reductions$umap@cell.embeddings
# write.table(xpca_castrated_integrated, file = "Cell_Ranger_output/xpca_castrated_integrated.csv", sep = "\t")
# write.table(xumap_castrated_integrated, file = "Cell_Ranger_output/xumap_castrated_integrated.csv", sep = "\t")

sce <- SingleCellExperiment(assays = List(counts = t(X_castrated)))

# Gene filtering
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

# Normalization + calculate PCA
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)
# pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
# rd1 <- pca$x[,1:30]
# reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = xumap_castrated)
reducedDims(sce) <- SimpleList(integrated.PCA = xpca_castrated_integrated, integrated.UMAP = xumap_castrated_integrated, UMAP = xumap_castrated)

# annotation
colData(sce)$annotation <- meta_castrated$annot.transfer
colData(sce)$annot.id[colData(sce)$annotation == "Basal"] <- 1
colData(sce)$annot.id[colData(sce)$annotation == "Luminal"] <- 2
colData(sce)$annot.id[colData(sce)$annotation == "Mesenchymal"] <- 3

# run on the integrated space
# infer pt for pp and ppc separately
sce.pp <- sce[,as.vector(meta_castrated["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")]
sce.ppc <- sce[,as.vector(meta_castrated["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")]

sce.pp <- slingshot(sce.pp, clusterLabels = 'annotation', reducedDim = 'UMAP', start.clus = "Basal")
ss.pp <- colData(sce.pp)$slingshot
pt.pp <- ss.pp@assays@data@listData$pseudotime

sce.ppc <- slingshot(sce.ppc, clusterLabels = 'annotation', reducedDim = 'UMAP', start.clus = "Basal")
ss.ppc <- colData(sce.ppc)$slingshot
pt.ppc <- ss.ppc@assays@data@listData$pseudotime

meta_castrated[meta_castrated["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)","pt_slingshot"] <- pt.pp[,"Lineage1"]
meta_castrated[meta_castrated["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)","pt_slingshot"] <- pt.ppc[,"Lineage1"]

write.table(meta_castrated, file = "results_ti/castrated_basal_mesenchymal/epithelial_mesenchymal_pt.csv", sep = "\t", quote = FALSE)

# Downstream analysis
library(tradeSeq)
# fit negative binomial GAM
genes <- rownames(sce)
sce.pp <- fitGAM(sce.pp, genes = genes)
sce.ppc <- fitGAM(sce.ppc, genes = genes)

# test for dynamic expression
ATres.pp <- associationTest(sce.pp)
write.table(ATres.pp, "results_ti/castrated_basal_mesenchymal/de_ti_tradeseq_pp.csv", sep = "\t", quote = FALSE)

ATres.ppc <- associationTest(sce.ppc)
write.table(ATres.ppc, "results_ti/castrated_basal_mesenchymal/de_ti_tradeseq_ppc.csv", sep = "\t", quote = FALSE)

# select the top genes with the smallest pvalues
topgenes <- rownames(ATres.pp[order(ATres.pp$pvalue), ])[1:50]
# order the cells according to the pseudotime 
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$annot.id[pst.ord]

heatmap(log1p(data.matrix(heatdata)), Colv = NA, ColSideColors = brewer.pal(9,"Set1")[heatclus])
