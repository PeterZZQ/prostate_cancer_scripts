rm(list=ls())
gc()

setwd("/localscratch/ziqi/Prostate_Cancer")
library("Seurat")
library("anndata")
library("SeuratDisk")
library("patchwork")

prostate.reference.h5ad <- read_h5ad("reference_data/reference_raw.h5ad")
# select batch T00_intact_1
prostate.reference.h5ad <- prostate.reference.h5ad[prostate.reference.h5ad$obs$batchID == "T00_intact_1",]
prostate.reference.obj <- CreateSeuratObject(counts = t(prostate.reference.h5ad$X), meta.data = prostate.reference.h5ad$obs)
prostate.reference.obj <- NormalizeData(prostate.reference.obj)
prostate.reference.obj <- FindVariableFeatures(prostate.reference.obj, selection.method = "vst", nfeatures = 2000)
# IntegrateData
DefaultAssay(prostate.reference.obj) <- "RNA"
# Run the standard workflow for visualization and clustering
prostate.reference.obj <- ScaleData(prostate.reference.obj, verbose = FALSE)
prostate.reference.obj <- RunPCA(prostate.reference.obj, npcs = 30, verbose = FALSE)
prostate.reference.obj <- RunUMAP(prostate.reference.obj, reduction = "pca", dims = 1:30, verbose = FALSE)
p1 <- DimPlot(prostate.reference.obj, reduction = "umap", group.by = "IntType")
p2 <- DimPlot(prostate.reference.obj, reduction = "umap", group.by = "batchID")
p1 + p2

prostate.h5ad <- read_h5ad("Cell_Ranger_output/adata_merge.h5ad")
prostate.obj <- CreateSeuratObject(counts = t(prostate.h5ad$X), meta.data = prostate.h5ad$obs)
# split the dataset into a list of seurat objects by batches
prostate.obj.list <- SplitObject(prostate.obj, split.by = "sample")

# normalize and identify variable features for each dataset independently
prostate.obj.list <- lapply(X = prostate.obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x.anchors <- FindTransferAnchors(reference = prostate.reference.obj, query = x, dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = x.anchors, refdata = prostate.reference.obj$IntType, dims = 1:30)
  x <- AddMetaData(x, metadata = predictions)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = prostate.obj.list)
prostate.anchors <- FindIntegrationAnchors(object.list = prostate.obj.list, anchor.features = features)
prostate.combined <- IntegrateData(anchorset = prostate.anchors)
DefaultAssay(prostate.combined) <- "integrated"

prostate.combined <- ScaleData(prostate.combined, verbose = FALSE)
prostate.combined <- RunPCA(prostate.combined, npcs = 30, verbose = FALSE)
prostate.combined <- RunUMAP(prostate.combined, reduction = "pca", dims = 1:30, n.neighbors = 50L, min.dist = 0.3)

# ------------------export the predicted and cluster label------------------
SaveH5Seurat(prostate.combined, filename = "Cell_Ranger_output/seurat_transferred_label.h5Seurat", overwrite = TRUE)
prostate.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_transferred_label.h5Seurat")

prostate.combined2 <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
meta.data <- prostate.combined@meta.data
meta.data$seurat_clusters <- prostate.combined2@meta.data$seurat_clusters
prostate.h5ad$obs$predicted_id <- meta.data[rownames(prostate.h5ad$obs), "predicted.id"] 
prostate.h5ad$obs$seurat_clusters <- meta.data[rownames(prostate.h5ad$obs), "seurat_clusters"]
write_h5ad(prostate.h5ad, filename = "Cell_Ranger_output/adata_seurat.h5ad")
# --------------------------------------------------------------------------

# Visualization
p1 <- DimPlot(prostate.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(prostate.combined, reduction = "umap", group.by = "predicted.id")
p1 + p2
ggsave(filename = "integration/transfer_label.pdf", device = "pdf", plot = p1 + p2, width = 12, height = 5)
ggsave(filename = "integration/transfer_label.png", device = "png", plot = p1 + p2, width = 12, height = 5, dpi = 150)


prostate.m1416_pp18 <- prostate.combined[,prostate.combined$sample == "M1416-PP18"]
prostate.m1417_pp12 <- prostate.combined[,prostate.combined$sample == "M1417-PP12"]
prostate.m1436_ppc12 <- prostate.combined[,prostate.combined$sample == "M1436-PPC12"]
prostate.m1437_ppc18 <- prostate.combined[,prostate.combined$sample == "M1437-PPC18"]

p1 <- DimPlot(prostate.m1416_pp18, reduction = "umap", group.by = "predicted.id")
p2 <- DimPlot(prostate.m1417_pp12, reduction = "umap", group.by = "predicted.id")
p3 <- DimPlot(prostate.m1436_ppc12, reduction = "umap", group.by = "predicted.id")
p4 <- DimPlot(prostate.m1437_ppc18, reduction = "umap", group.by = "predicted.id")
p1 + p2 + p3 + p4
ggsave(filename = "integration/transfer_label_separate.pdf", device = "pdf", plot = p1 + p2 + p3 + p4, width = 12, height = 8)
ggsave(filename = "integration/transfer_label_separate.png", device = "png", plot = p1 + p2 + p3 + p4, width = 12, height = 8, dpi = 150)

