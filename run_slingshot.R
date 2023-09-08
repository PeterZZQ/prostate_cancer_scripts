rm(list = ls())
gc()
library(anndata)
library(SingleCellExperiment)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(slingshot)
library(umap)

setwd("/localscratch/ziqi/Prostate_Cancer")
ad <- read_h5ad("Cell_Ranger_output/adata_ti.h5ad")
ad_m1417_pp12 <- ad[ad$obs$sample == "M1417-PP12",]
ad_m1416_pp18 <- ad[ad$obs$sample == "M1416-PP18",]
ad_m1436_ppc12 <- ad[ad$obs$sample == "M1436-PPC12",]
ad_m1437_ppc18 <- ad[ad$obs$sample == "M1437-PPC18",]

X_pca <- ad$obsm$X_pca
counts <- t(ad$X)
sce <- SingleCellExperiment(assays = List(counts = counts))
pca <- prcomp(t(counts), scale. = FALSE)
X_pca <- pca$x[,1:30]
reducedDims(sce) <- SimpleList(PCA = X_pca)
colData(sce)$annotation <- ad$obs$annotation
colData(sce)$leiden <- ad$obs$leiden
sce <- slingshot(sce, clusterLabels = 'leiden', reducedDim = 'PCA', start.clus = "2")
ss <- colData(sce)$slingshot
pt <- ss@assays@data@listData$pseudotime
ad$obs[,"pseudotime_joint"] <- pt[,"Lineage1"]


X_pca <- ad_m1417_pp12$obsm$X_pca
counts <- t(ad_m1417_pp12$X)
sce <- SingleCellExperiment(assays = List(counts = counts))
pca <- prcomp(t(counts), scale. = FALSE)
X_pca <- pca$x[,1:30]
reducedDims(sce) <- SimpleList(PCA = X_pca)
colData(sce)$annotation <- ad_m1417_pp12$obs$annotation
colData(sce)$leiden <- ad_m1417_pp12$obs$leiden
sce <- slingshot(sce, clusterLabels = 'leiden', reducedDim = 'PCA', start.clus = "2")
ss <- colData(sce)$slingshot
pt_m1417_pp12 <- ss@assays@data@listData$pseudotime


X_pca <- ad_m1416_pp18$obsm$X_pca
counts <- t(ad_m1416_pp18$X)
sce <- SingleCellExperiment(assays = List(counts = counts))
pca <- prcomp(t(counts), scale. = FALSE)
X_pca <- pca$x[,1:30]
reducedDims(sce) <- SimpleList(PCA = X_pca)
colData(sce)$annotation <- ad_m1416_pp18$obs$annotation
colData(sce)$leiden <- ad_m1416_pp18$obs$leiden
sce <- slingshot(sce, clusterLabels = 'leiden', reducedDim = 'PCA', start.clus = "2")
ss <- colData(sce)$slingshot
pt_m1416_pp18 <- ss@assays@data@listData$pseudotime


X_pca <- ad_m1436_ppc12$obsm$X_pca
counts <- t(ad_m1436_ppc12$X)
sce <- SingleCellExperiment(assays = List(counts = counts))
pca <- prcomp(t(counts), scale. = FALSE)
X_pca <- pca$x[,1:30]
reducedDims(sce) <- SimpleList(PCA = X_pca)
colData(sce)$annotation <- ad_m1436_ppc12$obs$annotation
colData(sce)$leiden <- ad_m1436_ppc12$obs$leiden
sce <- slingshot(sce, clusterLabels = 'leiden', reducedDim = 'PCA', start.clus = "2")
ss <- colData(sce)$slingshot
pt_m1436_ppc12 <- ss@assays@data@listData$pseudotime


X_pca <- ad_m1437_ppc18$obsm$X_pca
counts <- t(ad_m1437_ppc18$X)
sce <- SingleCellExperiment(assays = List(counts = counts))
pca <- prcomp(t(counts), scale. = FALSE)
X_pca <- pca$x[,1:30]
reducedDims(sce) <- SimpleList(PCA = X_pca)
colData(sce)$annotation <- ad_m1437_ppc18$obs$annotation
colData(sce)$leiden <- ad_m1437_ppc18$obs$leiden
sce <- slingshot(sce, clusterLabels = 'leiden', reducedDim = 'PCA', start.clus = "2")
ss <- colData(sce)$slingshot
pt_m1437_ppc18 <- ss@assays@data@listData$pseudotime

ad$obs$pseudotime <- 0
ad$obs[ad$obs$sample == "M1417-PP12", "pseudotime"] <- pt_m1417_pp12
ad$obs[ad$obs$sample == "M1416-PP18", "pseudotime"] <- pt_m1416_pp18
ad$obs[ad$obs$sample == "M1436-PPC12", "pseudotime"] <- pt_m1436_ppc12
ad$obs[ad$obs$sample == "M1437-PPC18", "pseudotime"] <- pt_m1437_ppc18

write_h5ad(ad, filename = "Cell_Ranger_output/adata_ti.h5ad")
